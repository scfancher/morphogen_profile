clear all
close all
clc

NVS = [35 40 90 100];
LVS = linspace(.01,25,1000);
PVS = logspace(-3,1,1000);

CN = zeros(length(NVS),length(LVS),length(PVS));

for i = 1:length(NVS)
    disp(i)
    for j = 1:length(LVS)
        sNDR = 2*exp(2*(1:NVS(i)-1)/LVS(j))*(1-exp(-1/LVS(j))*(2/LVS(j)+sinh(2/LVS(j)))/...
            (4*sinh(1/LVS(j))))/(4*exp(-2/LVS(j))*(sinh(1/LVS(j)))^3);
        for k = 1:length(PVS)
            gam = gammavals(NVS(i),LVS(j),PVS(k));
            cNDR = 2*gam(NVS(i)+1:end-1)*sum(gam)./((gam(NVS(i)+1:end-1)-gam(NVS(i)+2:end)).^2);
            rho = sNDR./cNDR;
            CN(i,j,k) = length(find(rho(:)<1 & isnan(rho(:))==0));
        end
    end
    for l = 0:NVS(i)-1
        check = 1;
        for j = 1:length(LVS)
           check = check*(min(CN(i,j,:))~=l);
        end
        if check==1
            fprintf('NO CN=%d for N=%d\n',l,NVS(i));
        end
    end
end

save('CNdata.mat','NVS','LVS','PVS','CN');

figure(1)
hold on
for j = 1:length(LVS)
    CNP = zeros(size(PVS));
    for k = 1:length(PVS)
        CNP(k) = CN(4,j,k);
    end
    plot(PVS,CNP)
end
hold off
set(gca,'Xscale','log')

figure(2)
hold on
for i = 1:length(NVS)
    CNL = zeros(size(LVS));
    for j = 1:length(LVS)
        CNL(j) = min(CN(i,j,:));
    end
    plot(LVS,CNL)
end
hold off