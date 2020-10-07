clear all
close all
clc

N = 50;
lambda = 10;

phivals = logspace(-3,1);
ndrvals = zeros(N-1,length(phivals));
ndrleg = cell(1,N-1);
for i = 1:N-1
    for j = 1:length(phivals)
        ndrvals(i,j) = cytoNDR(i,N,lambda,phivals(j));
    end
    ndrleg{i} = sprintf('j=%d',i);
end

phic = zeros(1,N-1);
for i = 1:N-1
    phic(i) = fminsearch(@(p) cytoNDR(i,N,lambda,p),0.5);
end

figure(1)
loglog(phivals,ndrvals)
%legend(ndrleg)
xlabel('\phi')
ylabel('NDR*\beta T')

figure(2)
semilogy(1:N-1,phic)
xlabel('j')
ylabel('\phi_c')