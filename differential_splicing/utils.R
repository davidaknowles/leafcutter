

# logistic (sigmoid) function
logistic=function(g) 1/(1+exp(-g))

# get inverse of PSD matrix using eigendecomposition
# make robust by thresholding all eigenvalues to be > eigenThreshold
robustSolve=function(si, eigenThreshold=0.01){
  svdSigma=eigen(si)
  svdSigma$values=pmax(svdSigma$values,eigenThreshold)
  svdSigma$vectors %*% diag(1/svdSigma$values) %*% t(svdSigma$vectors)
}