function ndr = cytoNDR(j,N,lambda,phi)
    gam = gammavals(N,lambda,phi);
    ndr = 2*gam(N+j)*sum(gam)/((gam(N+j)-gam(N+j+1))^2);
end