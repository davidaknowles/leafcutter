import torch

import pyro
import pyro.distributions as dist
from torch.distributions import constraints

import numpy as np
import scipy.stats
import time
from sklearn import linear_model

from leafcutter.differential_splicing.optim import fit_with_lbfgs

from dataclasses import dataclass

@dataclass
class LeafcutterFit:
    """
    Stores the result of fitting the LeafCutter differential splicing model to one cluster. 
    """
    beta: torch.Tensor
    conc: torch.Tensor
    loss: float
    exit_status: str

class LeafCutterModel(pyro.nn.PyroModule):
    """
    The LeafCutter differential splicing model for one cluster. 
    """
    
    def __init__(self, P, J, eps = 1e-8, gamma_shape = 1.0001, gamma_rate = 1e-4, beta_scale = np.inf, multiconc = True): 
        """
        Initialize the LeafCutterModel.

        Args:
            P (int): Number of covariates.
            J (int): Number of junctions.
            eps (float): A small constant to prevent numerical issues.
            gamma_shape (float): Shape parameter for the gamma distribution used for the concentration prior. 
            gamma_rate (float): Rate parameter for the gamma distribution used for the concentration prior. 
            beta_scale (float): prior std for coefficients. 
            multiconc (bool): Indicates whether to use a separate concentration parameter for each junction.

        """
        super().__init__()
        self.P = P
        self.J = J
        self.eps = eps
        self.multinomial = gamma_shape is None
        if not self.multinomial: 
            conc_prior = dist.Gamma(gamma_shape,gamma_rate)
            self.conc_prior = conc_prior.expand([self.J]).to_event(1) \
                if multiconc else conc_prior
        self.beta_scale = beta_scale
       
    def forward(self, x, y): 
        """
        Run the generative process. 

        Args:
            x (torch.Tensor): Input data representing covariates.
            y (torch.Tensor): Observed data, i.e. junction counts. 
        """
        
        with pyro.plate("covariates", self.P): # beta is P (covariates) x J (junctions)
            beta_dist = dist.Normal(0.0,self.beta_scale) if np.isfinite(self.beta_scale) else dist.ImproperUniform(constraints.real, (), ())
            b_param = pyro.sample("beta", beta_dist.expand([self.P, self.J]).to_event(1)) 

        logits = x @ b_param
        
        if self.multinomial: 
            logits_norm = logits - logits.logsumexp(1, keepdims = True)
            pyro.factor("multinomial", (y * logits_norm).sum())
        else:
            conc_param = pyro.sample("conc", self.conc_prior)
            a = logits.softmax(1) * conc_param + self.eps
            sum_a = a.sum(1)
            
            # manual implementation of the likelihood is faster since it avoids calculating normalization constants
            pyro.factor("dm", (sum_a.lgamma() + (a+y).lgamma().sum(1) - (sum_a + y.sum(1)).lgamma() - a.lgamma().sum(1)).sum())
            
            #y_dist = dist.DirichletMultinomial(a, total_count = y.sum(1))
            #with pyro.plate("data", y.shape[0]):
            #    return pyro.sample("obs", y_dist, obs=y)

class BaseGuide():
    """
    Has common functionality for the different possible guides, shouldn't be used directly. 
    """
    
    def __init__(self, init_beta, multinomial = False, init_conc = None, multiconc = True, conc_max = 300.):
        self.multiconc = multiconc
        self.conc_max = conc_max
        (P,J) = init_beta.shape
        torch_types = { "device" : init_beta.device, "dtype" : init_beta.dtype }
        self.multinomial = multinomial
        if not multinomial: 
            default_init_conc = torch.full([J],10.,**torch_types) if multiconc else torch.tensor(10.,**torch_types)
            self.init_conc = default_init_conc if (init_conc is None) else init_conc

    @property
    def conc(self):
        if self.multinomial: return torch.inf
        return pyro.param("conc_loc", lambda: self.init_conc.clone().detach(), constraint=constraints.interval(0., self.conc_max))
    
    @property
    def beta(self):
        raise NotImplementedError("Subclasses should implement this!")
    
    def __call__(self, x, y): # ok i think this acutally incorporates the constraint
        if not self.multinomial: 
            conc_param = pyro.param("conc_loc", lambda: self.init_conc.clone().detach(), constraint=constraints.interval(0., self.conc_max))
            pyro.sample("conc", dist.Delta(self.conc, event_dim = 1 if self.multiconc else 0)) # do we need to adjust for the constraint manually here? 
        beta_param = self.beta
        with pyro.plate("covariates", beta_param.shape[0]):
            pyro.sample("beta", dist.Delta(beta_param, event_dim = 1))

        
class SimpleGuide(BaseGuide):
    """
    Doesn't deal with the extra degree of freedom in beta. 
    """
    
    def __init__(self, init_beta, **kwargs):
        super().__init__(init_beta, **kwargs)
        self.init_beta = init_beta

    @property
    def beta(self):
        return pyro.param("beta_loc", lambda: self.init_beta.clone().detach())

class CleverGuide(BaseGuide):
    """
    This is the parametrization that the original LeafCutter R/Stan package used. 
    """
    
    def __init__(self, init_beta, **kwargs):
        super().__init__(init_beta, **kwargs)
        (P,J) = init_beta.shape

        up = init_beta[ [range(P),init_beta.abs().max(1).indices] ] / (1.-1./J)
        down = init_beta[ [range(P),(-init_beta * up.sign()[:,None]).max(1).indices] ] / (1./J)
        beta_scale = up - down
        beta_raw = init_beta / beta_scale[:,None] + 1./J
        beta_raw[beta_scale == 0.,:] = 1./J # can this learn?
        assert(not beta_raw.isnan().any().item())
        beta_raw = beta_raw / beta_raw.sum(1, keepdim = True)
        assert(not beta_raw.isnan().any().item())
        
        self.init_beta_raw = beta_raw
        self.init_beta_scale = beta_scale

        #assert( (init_beta - (beta_raw - 1./J) * beta_scale[:,None]).abs().max().item() < 1e-5 )
        #assert( ((beta_raw.sum(1) - 1.).abs().mean() < 1e-6).item() )
        #beta_reconstructed = beta_scale[:,None] * (beta_raw - 1./J)
        #assert( ((beta_reconstructed - init_beta).abs().mean() < 1e-5).item() )

    @property
    def beta(self):
        beta_raw_param = pyro.param("beta_raw", lambda: self.init_beta_raw.clone().detach(), constraint=constraints.simplex)
        beta_scale_param = pyro.param("beta_scale", lambda: self.init_beta_scale.clone().detach())
        return beta_scale_param[:,None] * (beta_raw_param - 1./beta_raw_param.shape[1]) 

class DamCleverGuide(BaseGuide):
    """
    This is the more recently proposed approach in the Stan docs: https://mc-stan.org/docs/stan-users-guide/parameterizing-centered-vectors.html under `QR decomposition`. Originally proposed by Aaron J Goodman. 
    """

    def __init__(self, init_beta, **kwargs):
        super().__init__(init_beta, **kwargs)
        J = init_beta.shape[1]
        
        A = torch.eye(J, dtype = init_beta.dtype, device = init_beta.device)
        A[-1,:-1] = -1.
        A[-1,-1] = 0. 
        self.A_qr = torch.linalg.qr(A).Q[:,:-1] # [J x (J-1)]
        # v = torch.eye(J-1) / (1.-1./J)
        # A_qr @ v @ A_qr.t() # gives correct marginals! 
        # beta_trans = A_qr @ (torch.randn(J-1) / np.sqrt(1.-1./J)) # sums to 0! 
        self.init_beta_raw = torch.linalg.solve(self.A_qr.t() @ self.A_qr, self.A_qr.t() @ init_beta.t()).t()

    @property
    def beta(self):
        beta_raw_param = pyro.param("beta_raw", lambda: self.init_beta_raw.clone().detach()) # note this transform does not require a Jacobian since it is constant/linear
        return beta_raw_param @ self.A_qr.t()

def brr_initialization(x, y): 
    """
    Try to get a good initialization using Bayesian ridge regression per junction. 
    """
    y_norm = torch.log( (y+1) / (y+1).sum(1, keepdim = True) ).cpu().numpy()
    x_np = x.cpu().numpy()
    
    N,P = x.shape
    J = y_norm.shape[1]
    beta_mm = np.zeros([P,J])
    reg = linear_model.BayesianRidge()
    for j in np.arange(y.shape[1]): 
        reg.fit(x_np, y_norm[:,j])
        beta_mm[:,j] = reg.coef_
    beta_mm = torch.tensor(beta_mm, dtype = x.dtype, device = x.device)
    
    return beta_mm - beta_mm.mean(1, keepdim = True)

def rr_initialization(x, y, regularizer = 0.001): 
    """
    Try to get a good initialization for beta by moment matching. V slightly slower than BRR overall.
    """
    y_norm = torch.log( (y+1) / (y+1).sum(1, keepdim = True) )
    # get estimate by moment matching
    I = torch.eye(x.shape[1], dtype = x.dtype, device = x.device)
    beta_mm = torch.linalg.solve( x.t() @ x + regularizer * I, x.t() @ y_norm )
    
    return beta_mm - beta_mm.mean(1, keepdim = True)

def fit_multinomial_glm(x, y, beta_init = None, fitter = fit_with_lbfgs, guide_type = DamCleverGuide): 
    """
    Try to get a good initialization for beta by moment matching. V slightly slower than BRR overall.
    """
    [N,P]=x.shape
    J = y.shape[1]
    pyro.clear_param_store()
    
    if beta_init is None: 
        beta_init = torch.zeros(P, J, device = x.device, dtype = x.dtype)

    multinomial_model = LeafCutterModel(P, J, None, gamma_shape = None)
    guide = guide_type(beta_init, multinomial = True)
    losses = fitter(multinomial_model, guide, x, y)
    return guide.beta

def fit_dm_glm(x, y, beta_init, conc_max = 3000., concShape=1.0001, concRate=1e-4, multiconc = True, fitter = fit_with_lbfgs, guide_type = DamCleverGuide, eps = 1.0e-8): 
    """
    Fits a Dirichlet Multinomial generalized linear model.

    Args:
        x (torch.Tensor): The input data with shape (N, P) where N is the number of observations and P is the number of predictors.
        y (torch.Tensor): The target data with shape (N, J) where J is the number of splice junctions.
        beta_init (torch.Tensor): The initial values for the coefficients, shape=(P,J)
        conc_max (float, optional): Maximum concentration parameter value. Defaults to 300.
        concShape (float, optional): Shape parameter for concentration prior. Defaults to 1.0001.
        concRate (float, optional): Rate parameter for concentration prior. Defaults to 1e-4.
        multiconc (bool, optional): Whether to use multiple concentration parameters. Defaults to True.
        fitter (function, optional): The fitting function to use. Defaults to fit_with_lbfgs.
        guide_type (class, optional): The type of guide to use. Defaults to DamCleverGuide.
        eps (float): small pseudocount added to DM parameter for numerical stability. 

    Returns:
        LeafcutterFit: An object containing fitted model parameters, losses, and exit status.
    """
    pyro.clear_param_store()

    (N,P) = x.shape
    J = y.shape[1]
    
    model = LeafCutterModel(P, J, gamma_shape = concShape, gamma_rate = concRate, multiconc = multiconc, eps = eps)
    guide = guide_type(beta_init, multiconc = multiconc, conc_max = conc_max)
    losses, exit_status = fitter(model, guide, x, y)

    return LeafcutterFit(
            beta = guide.beta.clone().detach(), 
            conc = guide.conc.clone().detach(),
            loss = losses[-1],
            exit_status = exit_status
        ) 

def simple_simulation(N, P, J, total_count = 100, conc = 10.):
    """
    Very simple simulation of data for one cluster. 
    
    Args:
        N (int): Number of samples.
        P (int): Number of covariates.
        J (int): Number of junctions.
        total_count (int): Total count for the cluster (default: 100).
        conc (float): Concentration parameter for the Dirichlet-Multinomial distribution (default: 10.0).

    Returns:
        tuple: A tuple containing the following elements:
            - x (torch.Tensor): Simulated covariate data.
            - y (torch.Tensor): Simulated observed data.
            - true_beta_norm (torch.Tensor): True beta parameters after normalization.
            - g (torch.Tensor): Probabilities computed using the softmax function.

    """
    x = torch.randn((N,P))
    x[:,0] = 1. # intercept
    b = torch.randn((P,J))
    xb = x @ b
    g = torch.softmax( x @ b, 1 ) 
    true_beta_norm = b - b.mean(1, keepdim = True)
    dm = dist.DirichletMultinomial(g * conc, total_count = total_count)
    y = dm.sample()
    return(x,y,true_beta_norm,g)


def dirichlet_multinomial_anova(x_full, x_null, y, init = "brr", **kwargs): 
    """
    Perform Dirichlet-Multinomial ANOVA analysis.

    Args:
        x_full (torch.Tensor): Full (i.e. alternative hypothesis) covariate data.
        x_null (torch.Tensor): Null (i.e. null hypothesis) covariate data.
        y (torch.Tensor): Observed junction counts
        init (str): initialization strategy. One of "brr" (Bayesian ridge regression), "rr" (ridge regression), "mult" (multinomial logistic regression) or "0" (set to 0). 
        concShape (float): Shape parameter for concentration priors (default: 1.0001).
        concRate (float): Rate parameter for concentration priors (default: 1e-4).
        multiconc (bool): Indicates whether to use separate concentration parameters for each junction (default: False).
        fitter (function): A function used to fit the model (default: fit_with_lbfgs).
        guide_type: one of BasicGuide, CleverGuide or DamCleverGuide

    Returns:
        tuple: A tuple containing the following elements:
            - loglr (float): Log-likelihood ratio statistic.
            - df (int): Degrees of freedom.
            - lrtp (float): Likelihood ratio test p-value.
            - null_fit (object): Result of fitting the null model.
            - full_fit (object): Result of fitting the full model.
            - refit_null_flag (bool): Flag indicating if the null model was refitted based on the full model.
    """
    
    (N,P_full) = x_full.shape
    (N,P_null) = x_null.shape
    J = y.shape[1]
    
    torch_types = { "device" : y.device, "dtype" : y.dtype }
    
    def get_init(x,y): 
        if init == "brr":
            return brr_initialization(x,y)
        elif init == "rr": 
            return rr_initialization(x,y)
        elif init == "mult": # doesn't speed things up
            mult_kwargs = { k:v for k,v in kwargs.items() if k in ["fit_with_lbfgs", "guide_type"]}
            return fit_multinomial_glm(x_null, y, **mult_kwargs) 
        elif init == "0": 
            return torch.zeros(x.shape[1], J, **torch_types)
        else: 
            raise Exception(f"Unknown initialization strategy {init}")
   
    # fit null model
    start_time = time.time()
    beta_init_null = get_init(x_null,y)
    null_fit = fit_dm_glm( x_null, y, beta_init_null, **kwargs )
    elapsed_time = time.time() - start_time
    #print(f"Elapsed time: {elapsed_time} seconds, {null_fit.loss}")

    # fit full model, initialized at null model solution
    beta_init_full = torch.cat( (null_fit.beta, torch.zeros((P_full - P_null, J), **torch_types)) )
    #beta_init_full -= beta_init_full.mean(1, keepdim = True) # shouldn't be necessary with the clever guides
    full_fit = fit_dm_glm( x_full, y, beta_init_full, **kwargs )
    
    # fit full model, initialized "smartly", and check if this gives a better fit
    beta_init_full = get_init(x_full,y)
    full_fit_smart = fit_dm_glm( x_full, y, beta_init_full, **kwargs )
    if full_fit_smart.loss < full_fit.loss: # this initialization did better
        full_fit = full_fit_smart
    
    df=(P_full - P_null)*(J-1)
    refit_null_flag=False
    loglr = null_fit.loss - full_fit.loss
    lrtp = scipy.stats.chi2(df).sf(2.*loglr) # sf = 1-cdf
    
    if lrtp < 0.001: # if result looks highly significant, check if we could improve fit initializing null based on full
        beta_init_null = full_fit.beta[:P_null,:]
        refit_null = fit_dm_glm(x_null, y, beta_init_null, **kwargs )
        if refit_null.loss < refit_null.loss: # if new fit is better
            null_fit = refit_null
            refit_null_flag = True
            loglr = null_fit.loss-full_fit.loss
            lrtp = scipy.stats.chi2(df = df).sf(2.*loglr)
    
    return loglr, df, lrtp, null_fit, full_fit, refit_null_flag
