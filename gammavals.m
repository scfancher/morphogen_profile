function gam = gammavals(N,lambda,phi)
    kappa = sign(phi)*(exp(abs(phi))-1)*(abs(phi)-log(exp(abs(phi))-1))/lambda;
    gam1s = exp(-kappa*2*(1:N))*(1-exp(-phi))./(1-exp(-phi-kappa*2*(1:N)));
    gam = [gam1s(end:-1:1) gam1s];
end