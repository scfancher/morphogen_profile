function ndr = SDCNDR(j,lambda)
    ndr = 2*exp(2*j/lambda)*(1-exp(-1/lambda)*(2/lambda+sinh(2/lambda))/...
        (4*sinh(1/lambda)))/(4*exp(-2/lambda)*(sinh(1/lambda))^3);
end