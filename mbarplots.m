close all
clear all
clc

N = 100;
lambda = 2;
phivals = logspace(-3,1,10);
legents = cell(1,length(phivals));

gamvals = zeros(length(phivals),2*N+1);
for i = 1:length(phivals)
    gam = gammavals(N,lambda,phivals(i));
    gamvals(i,:) = [gam(1:N) 1 gam(N+1:end)]/sum(gam);
    legents{i} = sprintf('phi=%3.3e',phivals(i));
end

figure(1)
plot(1:N,gamvals(:,N+2:end))
legend(legents)
set(gca,'YScale','log')
%set(gca,'XScale','log')