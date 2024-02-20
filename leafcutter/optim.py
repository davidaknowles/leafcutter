import torch

import pyro
from torch.distributions import constraints
from pyro.infer import SVI, Trace_ELBO

import numpy as np

import time

#from torch.optim.lbfgs import _strong_wolfe
from torch.optim.lbfgs import _cubic_interpolate

def _strong_wolfe(obj_func,
                  x,
                  t,
                  d,
                  f,
                  g,
                  gtd,
                  c1=1e-4,
                  c2=0.9,
                  tolerance_change=1e-9,
                  max_ls=25):
    # ported from https://github.com/torch/optim/blob/master/lswolfe.lua
    d_norm = d.abs().max()
    g = g.clone(memory_format=torch.contiguous_format)
    # evaluate objective and gradient using initial step
    f_new, g_new = obj_func(x, t, d)
    ls_func_evals = 1
    gtd_new = g_new.dot(d)

    # bracket an interval containing a point satisfying the Wolfe criteria
    t_prev, f_prev, g_prev, gtd_prev = 0, f, g, gtd
    done = False
    ls_iter = 0
    while ls_iter < max_ls:
        # check conditions
        if f_new > (f + c1 * t * gtd) or (ls_iter > 1 and f_new >= f_prev):
            bracket = [t_prev, t]
            bracket_f = [f_prev, f_new]
            bracket_g = [g_prev, g_new.clone(memory_format=torch.contiguous_format)]
            bracket_gtd = [gtd_prev, gtd_new]
            break

        if abs(gtd_new) <= -c2 * gtd:
            bracket = [t]
            bracket_f = [f_new]
            bracket_g = [g_new]
            done = True
            break

        if gtd_new >= 0:
            bracket = [t_prev, t]
            bracket_f = [f_prev, f_new]
            bracket_g = [g_prev, g_new.clone(memory_format=torch.contiguous_format)]
            bracket_gtd = [gtd_prev, gtd_new]
            break

        # interpolate
        min_step = t + 0.01 * (t - t_prev)
        max_step = t * 10
        tmp = t
        t = _cubic_interpolate(
            t_prev,
            f_prev,
            gtd_prev,
            t,
            f_new,
            gtd_new,
            bounds=(min_step, max_step))

        # next step
        t_prev = tmp
        f_prev = f_new
        g_prev = g_new.clone(memory_format=torch.contiguous_format)
        gtd_prev = gtd_new
        f_new, g_new = obj_func(x, t, d)
        ls_func_evals += 1
        gtd_new = g_new.dot(d)
        ls_iter += 1

    # reached max number of iterations?
    if ls_iter == max_ls:
        bracket = [0, t]
        bracket_f = [f, f_new]
        bracket_g = [g, g_new]

    # zoom phase: we now have a point satisfying the criteria, or
    # a bracket around it. We refine the bracket until we find the
    # exact point satisfying the criteria
    insuf_progress = False
    # find high and low points in bracket
    low_pos, high_pos = (0, 1) if bracket_f[0] <= bracket_f[-1] else (1, 0)
    while not done and ls_iter < max_ls:
        # line-search bracket is so small
        if abs(bracket[1] - bracket[0]) * d_norm < tolerance_change:
            break

        # compute new trial value
        t = _cubic_interpolate(bracket[0], bracket_f[0], bracket_gtd[0],
                               bracket[1], bracket_f[1], bracket_gtd[1])

        # test that we are making sufficient progress:
        # in case `t` is so close to boundary, we mark that we are making
        # insufficient progress, and if
        #   + we have made insufficient progress in the last step, or
        #   + `t` is at one of the boundary,
        # we will move `t` to a position which is `0.1 * len(bracket)`
        # away from the nearest boundary point.
        eps = 0.1 * (max(bracket) - min(bracket))
        if min(max(bracket) - t, t - min(bracket)) < eps:
            # interpolation close to boundary
            if insuf_progress or t >= max(bracket) or t <= min(bracket):
                # evaluate at 0.1 away from boundary
                if abs(t - max(bracket)) < abs(t - min(bracket)):
                    t = max(bracket) - eps
                else:
                    t = min(bracket) + eps
                insuf_progress = False
            else:
                insuf_progress = True
        else:
            insuf_progress = False

        # Evaluate new point
        f_new, g_new = obj_func(x, t, d)
        ls_func_evals += 1
        gtd_new = g_new.dot(d)
        ls_iter += 1

        if f_new > (f + c1 * t * gtd) or f_new >= bracket_f[low_pos]:
            # Armijo condition not satisfied or not lower than lowest point
            bracket[high_pos] = t
            bracket_f[high_pos] = f_new
            bracket_g[high_pos] = g_new.clone(memory_format=torch.contiguous_format)
            bracket_gtd[high_pos] = gtd_new
            low_pos, high_pos = (0, 1) if bracket_f[0] <= bracket_f[1] else (1, 0)
        else:
            if abs(gtd_new) <= -c2 * gtd:
                # Wolfe conditions satisfied
                done = True
            elif gtd_new * (bracket[high_pos] - bracket[low_pos]) >= 0:
                # old high becomes new low
                bracket[high_pos] = bracket[low_pos]
                bracket_f[high_pos] = bracket_f[low_pos]
                bracket_g[high_pos] = bracket_g[low_pos]
                bracket_gtd[high_pos] = bracket_gtd[low_pos]

            # new point becomes new low
            bracket[low_pos] = t
            bracket_f[low_pos] = f_new
            bracket_g[low_pos] = g_new.clone(memory_format=torch.contiguous_format)
            bracket_gtd[low_pos] = gtd_new

    # return stuff
    t = bracket[low_pos]
    f_new = bracket_f[low_pos]
    g_new = bracket_g[low_pos]
    return f_new, g_new, t, ls_func_evals

class MyLBFGS(torch.optim.LBFGS):
    """ Clone of torch.optim.LBFGS except that it 1) reports exit status and 2) returns a list of losses per iteration. """
    
    @torch.no_grad()
    def step(self, closure):
        """Perform a single optimization step.

        Args:
            closure (Callable): A closure that reevaluates the model
                and returns the loss.
        """
        assert len(self.param_groups) == 1

        # Make sure the closure is always called with grad enabled
        closure = torch.enable_grad()(closure)

        group = self.param_groups[0]
        lr = group['lr']
        max_iter = group['max_iter']
        max_eval = group['max_eval']
        tolerance_grad = group['tolerance_grad']
        tolerance_change = group['tolerance_change']
        line_search_fn = group['line_search_fn']
        history_size = group['history_size']

        # NOTE: LBFGS has only global state, but we register it as state for
        # the first param, because this helps with casting in load_state_dict
        state = self.state[self._params[0]]
        state.setdefault('func_evals', 0)
        state.setdefault('n_iter', 0)
        
        # evaluate initial f(x) and df/dx
        orig_loss = closure()
        loss = float(orig_loss)
        losses = [ loss ]
        
        current_evals = 1
        state['func_evals'] += 1

        flat_grad = self._gather_flat_grad()
        opt_cond = flat_grad.abs().max() <= tolerance_grad

        self.exit_status = "not set"
        
        # optimal condition
        if opt_cond:
            self.exit_status = "max|grad| < tolerance_grad at init"
            return orig_loss

        # tensors cached in state (for tracing)
        d = state.get('d')
        t = state.get('t')
        old_dirs = state.get('old_dirs')
        old_stps = state.get('old_stps')
        ro = state.get('ro')
        H_diag = state.get('H_diag')
        prev_flat_grad = state.get('prev_flat_grad')
        prev_loss = state.get('prev_loss')

        n_iter = 0
        # optimize for a max of max_iter iterations
        while n_iter < max_iter:
            # keep track of nb of iterations
            n_iter += 1
            state['n_iter'] += 1

            ############################################################
            # compute gradient descent direction
            ############################################################
            if state['n_iter'] == 1:
                d = flat_grad.neg()
                old_dirs = []
                old_stps = []
                ro = []
                H_diag = 1
            else:
                # do lbfgs update (update memory)
                y = flat_grad.sub(prev_flat_grad)
                s = d.mul(t)
                ys = y.dot(s)  # y*s
                if ys > 1e-10:
                    # updating memory
                    if len(old_dirs) == history_size:
                        # shift history by one (limited-memory)
                        old_dirs.pop(0)
                        old_stps.pop(0)
                        ro.pop(0)

                    # store new direction/step
                    old_dirs.append(y)
                    old_stps.append(s)
                    ro.append(1. / ys)

                    # update scale of initial Hessian approximation
                    H_diag = ys / y.dot(y)  # (y*y)

                # compute the approximate (L-BFGS) inverse Hessian
                # multiplied by the gradient
                num_old = len(old_dirs)

                if 'al' not in state:
                    state['al'] = [None] * history_size
                al = state['al']

                # iteration in L-BFGS loop collapsed to use just one buffer
                q = flat_grad.neg()
                for i in range(num_old - 1, -1, -1):
                    al[i] = old_stps[i].dot(q) * ro[i]
                    q.add_(old_dirs[i], alpha=-al[i])

                # multiply by initial Hessian
                # r/d is the final direction
                d = r = torch.mul(q, H_diag)
                for i in range(num_old):
                    be_i = old_dirs[i].dot(r) * ro[i]
                    r.add_(old_stps[i], alpha=al[i] - be_i)

            if prev_flat_grad is None:
                prev_flat_grad = flat_grad.clone(memory_format=torch.contiguous_format)
            else:
                prev_flat_grad.copy_(flat_grad)
            prev_loss = loss

            ############################################################
            # compute step length
            ############################################################
            # reset initial guess for step size
            if state['n_iter'] == 1:
                t = min(1., 1. / flat_grad.abs().sum()) * lr
            else:
                t = lr

            # directional derivative
            gtd = flat_grad.dot(d)  # g * d

            # directional derivative is below tolerance
            if gtd > -tolerance_change:
                self.exit_status = "Converged (directional derivative < tolerance_change)"
                break

            # optional line search: user function
            ls_func_evals = 0
            if line_search_fn is not None:
                # perform line search, using user function
                if line_search_fn != "strong_wolfe":
                    raise RuntimeError("only 'strong_wolfe' is supported")
                else:
                    x_init = self._clone_param()

                    def obj_func(x, t, d):
                        return self._directional_evaluate(closure, x, t, d)

                    loss, flat_grad, t, ls_func_evals = _strong_wolfe(
                        obj_func, x_init, t, d, loss, flat_grad, gtd)
                self._add_grad(t, d)
                opt_cond = flat_grad.abs().max() <= tolerance_grad
            else:
                # no line search, simply move with fixed-step
                self._add_grad(t, d)
                if n_iter != max_iter:
                    # re-evaluate function only if not in last iteration
                    # the reason we do this: in a stochastic setting,
                    # no use to re-evaluate that function here
                    with torch.enable_grad():
                        loss = float(closure())
                    flat_grad = self._gather_flat_grad()
                    opt_cond = flat_grad.abs().max() <= tolerance_grad
                    ls_func_evals = 1

            # update func eval
            current_evals += ls_func_evals
            state['func_evals'] += ls_func_evals
            
            losses.append( loss )

            ############################################################
            # check conditions
            ############################################################
            if n_iter == max_iter:
                self.exit_status = "Reached maximum # iterations (%i)" % n_iter
                break

            if current_evals >= max_eval:
                self.exit_status = "Reached maximum # evaluations (%i)" % current_evals
                break

            # optimal condition
            if opt_cond:
                self.exit_status = "Converged (line search optimality condition)" 
                break

            # lack of progress
            if d.mul(t).abs().max() <= tolerance_change:
                self.exit_status = "Converged (lack of progress)" 
                break

            if abs(loss - prev_loss) < tolerance_change:
                self.exit_status = "Converged (loss change below tolerance)" 
                break

        state['d'] = d
        state['t'] = t
        state['old_dirs'] = old_dirs
        state['old_stps'] = old_stps
        state['ro'] = ro
        state['H_diag'] = H_diag
        state['prev_flat_grad'] = prev_flat_grad
        state['prev_loss'] = prev_loss

        return losses # this implicitly return #iterations == len(losses)

def fit_with_SVI(model, guide, x, y, iterations = 1500, loss_tol = 1e-3, adam_opts = {"lr": 0.02}):
    """
    Fit the model using "SVI": but really there is no randomness so this should just be equivalent to stochastic gradient descent (SGD)
    """
    pyro.clear_param_store()
    optimizer = pyro.optim.Adam(adam_opts)
    svi = SVI(model, guide, optimizer, Trace_ELBO())
    losses = []
    old_loss = np.inf
    for i in range(iterations):
        #print(i,end="\r")
        loss = svi.step(x,y)
        losses.append(loss)
        if np.abs(loss - old_loss) < loss_tol: break
        old_loss = loss
    return losses

def manual_fit(model, guide, x, y, iterations = 1500, loss_tol = 1e-3, **adam_kwargs):
    """
    Direct implementation of SGD. 
    """
    pyro.clear_param_store()
    loss_fn = Trace_ELBO().differentiable_loss
    
    if not "lr" in adam_kwargs: adam_kwargs["lr"] = 0.02

    guide(x,y) # run once to instantiate params
    ps = pyro.get_param_store()
    params = [ p.unconstrained() for p in ps.values() ] # want to optimize in the constrained space 

    optimizer = torch.optim.Adam(params, **adam_kwargs)
    losses = []
    old_loss = np.inf
    for i in range(iterations):
        #print(i,end="\r")
        loss = loss_fn(model, guide, x, y)
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
        loss_item = loss.item()
        losses.append(loss_item)
        if np.abs(loss_item - old_loss) < loss_tol: break
        old_loss = loss_item
    return losses

def fit_with_lbfgs(model, guide, x, y, outer_iterations = 1, **lbfgs_kwargs):
    """
    Fit using LBFGS (which o.g. LeafCutter also used). 
    """
    my_defaults = { "max_iter" : 500, "line_search_fn" : "strong_wolfe" }
    for k,v in my_defaults.items(): 
        if not k in lbfgs_kwargs: 
            lbfgs_kwargs[k] = v
    
    pyro.clear_param_store()
    guide(x,y) # run once to instantiate params
    ps = pyro.get_param_store()
    params = [ p.unconstrained() for p in ps.values() ] # want to optimize in the constrained space 
    loss_fn = Trace_ELBO().differentiable_loss # JitTrace is actually slower it seems
    
    # defaults: lr=1, max_iter=20, max_eval=None, tolerance_grad=1e-07, tolerance_change=1e-09, history_size=100, line_search_fn=None
    optimizer = MyLBFGS(params, **lbfgs_kwargs) #lr=0.05, max_iter=inner_iterations, tolerance_grad=1e-4, history_size = 20)
    #optimizer = torch.optim.LBFGS(params, **lbfgs_kwargs)
    
    def closure(): # slightly awkward requirement of torch.optim.LBFGS
        optimizer.zero_grad()
        loss = loss_fn(model, guide, x, y)
        loss.backward()
        return loss
    losses = []
    for i in range(outer_iterations): # I think in principle we don't need this loop? 
        losses_here = optimizer.step(closure) # now returns a list of losses
        #losses += [ losses_here.item() ]
        losses += losses_here
    #losses.append(loss_fn(model, guide, x, y).item())
    return np.array(losses), optimizer.exit_status # "Converged" 
