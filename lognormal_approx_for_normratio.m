function [mu,sig]=lognormal_approx_for_normratio(mu1,sig1,mu2,sig2)
mu=fzero(@(x) diff_log_normratiopdf_of_exp(x,mu1,sig1,mu2,sig2),log(mu1/mu2));
[d,d2]=diff_log_normratiopdf_of_exp(mu,mu1,sig1,mu2,sig2);
sig=sqrt(-1./d2);

