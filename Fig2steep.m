close all
clear all
clc
load('CNdata');

mapN = 1000;
mymap = zeros(mapN,3);
for i = 1:floor(mapN/2)
    mymap(i,:) = [1 (i-1)/(floor(mapN/2)-1) (i-1)/(floor(mapN/2)-1)];
end
for i = floor(mapN/2)+1:mapN
    mymap(i,:) = [(mapN-i)/(mapN-floor(mapN/2)-1) (mapN-i)/(mapN-floor(mapN/2)-1) 1];
end

lamvals = [100 30 22.5 20.2 197 8 115 150 65 5.8];
errvals = [10 5 5 5.7 7 3 20 25 10 2.04];
avals = [2.8 10 2.8 1.32 10 1.32 10 10 10 1.32];
lamdata = lamvals./avals;
errdata = errvals./avals;
Ndata = [90 40 35 100 40 100 40 40 40 100];
mlabeldata = {'Bicoid' 'Cyclops' 'Dorsal' 'Dpp' 'Fgf8' 'Hh' 'Lefty1' 'Lefty2' 'Squint' 'Wg'};
source = [0 1 0 0 1 0 1 1 1 0];
mech = [1 0 0 0 1 -1 1 1 0 -1];

lamsort = zeros(size(lamvals));
ins = 0;
for i=-1:1
    in = find(mech(:)==i);
    [~,insort] = sort(lamvals(in));
    lamsort(ins+1:ins+length(in)) = in(insort);
    ins = ins+length(in);
end
lamhvals = lamdata(lamsort);
errsvals = errdata(lamsort);
Nvals = Ndata(lamsort);
sources = source(lamsort);
mechs = mech(lamsort);

f = figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

axA = subplot('Position',[.125 .55 .35 .4]);
lam1 = linspace(7,17,50);
N1 = 100;
j1 = 50;
rho1 = zeros(size(lam1));
for i = 1:length(lam1)
    phi1 = fminsearch(@(p) cytoNDR(j1,N1,lam1(i),p),0.1);
    rho1(i) = SDCNDR(j1,lam1(i))/cytoNDR(j1,N1,lam1(i),phi1);
end
semilogy(lam1,rho1,'k-')
hold on
plot([lam1(1) lam1(end)],[1 1],'k--')
hold off
xlim([lam1(1) lam1(end)])
xlabel('Profile length, $\hat{\lambda}=\lambda/a$','Interpreter','latex')
ylabel('$\rho_{j} = P_{DT}^{2}/P_{SDC}^{2}$','Interpreter','latex')
lamc = lam1(find((lam1-1).^2==min((lam1-1).^2),1));
phi1 = logspace(-2,0,50);
cNDR1 = zeros(size(phi1));
for i = 1:length(phi1)
    cNDR1(i) = cytoNDR(j1,N1,lamc,phi1(i));
end
ax = axA.Position;
axA_2 = axes('Position',[ax(1:2)+.5*ax(3:4) .425*ax(3:4)]);
loglog(phi1,1./cNDR1,'k-')
ylim([10^-7 10^-5.75])
xlabel('Shape parameter, $\phi$','Interpreter','latex','Fontsize',8)
ylabel('$P_{DT}^{2}/\beta T$','Interpreter','latex','Fontsize',8)
annotation('textbox',[0 .875 .1 .1],'String','A','LineStyle','none',...
    'Interpreter','latex','FontSize',16)

axC = subplot('Position',[.6 .55 .35 .05]);
nin100 = find(NVS(:)==100,1);
CNL = zeros(size(LVS));
for j = 1:length(LVS)
    CNL(j) = min(CN(nin100,j,:));
end
vertsC = zeros(200,2);
vertsC(1:2,:) = [0 0; 0 2];
facesC = zeros(99,4);
colsC = zeros(200,1);
colsC(1:2) = [0; 0];
for i = 2:99
    lambound = LVS(find(CNL(:)==i,1));
    vertsC(2*i-1:2*i,:) = [lambound 0; lambound 2];
    facesC(i-1,:) = [2*i-3 2*i-2 2*i 2*i-1];
    colsC(2*i-1:2*i) = [1; 1]*(i-2)/97;
    if i==50
        lamdiv = lambound;
        vertsC(end-1:end,:) = [2*lambound 0; 2*lambound 2];
        facesC(end,:) = [197 198 200 199];
        colsC(end-1:end) = [1; 1];
    end
end
vertsC(:,1) = vertsC(:,1)/lamdiv;
quants = quantile(vertsC(:,1),3);
hold on
patch('Faces', facesC, 'Vertices', vertsC, 'FaceVertexCData', colsC, 'FaceColor', 'interp', 'EdgeColor', 'none','FaceAlpha',.5)
plot([vertsC(3,1) vertsC(end-2,1)],[1 1],'k.')
plot([quants(1) quants(1)],[0 2],'k--')
plot([1 1],[0 2],'k-')
plot([quants(3) quants(3)],[0 2],'k--')
hold off
colormap(mymap)
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[])
xlabel('Profile length, $\hat{\lambda}/\hat{\lambda}_{50}$','Interpreter','latex')
axC_2 = axes('Position',axC.Position,'Color','none');
axC_2.XAxisLocation = 'top';
axC_2.XLim = axC.XLim;
axC_2.XTick = [vertsC(3,1),quants(1),1,quants(3),vertsC(end-2,1)];
axC_2.XTickLabel = {'0\%','25\%','50\%','75\%','100\%'};
axC_2.YTick = [];
xlabel('Fraction of cells more precise in SDC','Interpreter','latex')
box on
annotation('textbox',[.525 .585 .1 .1],'String','C','LineStyle','none',...
    'Interpreter','latex','FontSize',16)

axB = subplot('Position',[.6 .75 .35 .2]);
phicm = zeros(2,length(LVS));
for j = 1:length(LVS)
    phiset = PVS(find(CN(nin100,j,:)==CNL(j)));
    if isempty(phiset)==0
        phicm(1,j) = min(phiset);
        phicm(2,j) = max(phiset);
    end
end
%lamins = find(LVS(:)>vertsC(3,1)*lamdiv & LVS(:)<vertsC(end-2,1)*lamdiv);
lamins = 1:length(LVS);
semilogy(LVS(lamins)/lamdiv,phicm(1,lamins),'k-')
hold on
semilogy(LVS(lamins)/lamdiv,phicm(2,lamins),'k-')
fill([LVS(lamins) fliplr(LVS(lamins))]/lamdiv,[phicm(1,lamins) fliplr(phicm(2,lamins))],'r','FaceAlpha',0.5)
hold off
xlim([0 2])
xlabel('Profile length, $\hat{\lambda}/\hat{\lambda}_{50}$','Interpreter','latex')
ylim([10^-2 10^1])
ylabel('Optimal value, $\phi^{*}$','Interpreter','latex')
annotation('textbox',[.485 .875 .1 .1],'String','B','LineStyle','none',...
    'Interpreter','latex','FontSize',16)

axD = subplot('Position',[.125 .1 .35 .35]);
hold on
p1 = plot(NaN,NaN,'o','Markersize',6,'MarkerFaceColor','k','MarkerEdgeColor','k');
p2 = plot(NaN,NaN,'s','Markersize',6,'MarkerFaceColor','k','MarkerEdgeColor','k');
p3 = plot(NaN,NaN,'o','Markersize',6,'MarkerFaceColor','r','MarkerEdgeColor','k');
p4 = plot(NaN,NaN,'o','Markersize',6,'MarkerFaceColor','b','MarkerEdgeColor','k');
p5 = plot(NaN,NaN,'o','Markersize',6,'MarkerFaceColor','w','MarkerEdgeColor','k');
for i = 1:length(lamhvals)
    if sources(i) == 0
        marktype = 'o';
    elseif sources(i) == 1
        marktype = 's';
    end
    if mechs(i) == 1
        markcolor = 'b';
    elseif mechs(i) == 0
        markcolor = 'w';
    elseif mechs(i) == -1
        markcolor = 'r';
    end
    errorbar(lamvals(lamsort(i)),i,errvals(lamsort(i)),'horizontal',marktype,'Markersize',8,...
        'MarkerFaceColor',markcolor,'MarkerEdgeColor','k','Color','k','LineWidth',1,'CapSize',8)
end
hold off
legend([p1 p2 p3 p4 p5],{'{\it Drosophila}','Zebrafish','Evidence for DT','Evidence for SDC',...
    sprintf('Evidence for\nmultiple mechanisms')},'Interpreter','latex','Location','SouthEast','FontSize',6)
legend boxoff
xlabel('Profile length, $\lambda$ ($\mu$m)','Interpreter','latex')
ylim([0.5 length(lamhvals)+0.5])
yticks(1:length(lamhvals))
yticklabels(mlabeldata(lamsort))
box on
annotation('textbox',[0 .37 .1 .1],'String','D','LineStyle','none',...
    'Interpreter','latex','FontSize',16)

axE = subplot('Position',[.6 .1 .35 .35]);
verts = [];
faces = [];
col = [];
quants = zeros(length(lamhvals),5);
lamdivvec = zeros(size(lamhvals));
infill = 0;
for i = 1:length(lamhvals)
    nin100 = find(NVS(:)==Nvals(i),1);
    CNL = zeros(size(LVS));
    for j = 1:length(LVS)
        CNL(j) = min(CN(nin100,j,:));
    end
    vertsE = zeros(2*Nvals(i),2);
    vertsE(1:2,:) = [0 0; 0 2];
    facesE = zeros(Nvals(i)-1,4);
    colsE = zeros(2*Nvals(i),1);
    colsE(1:2) = [0; 0];
    for j = 2:Nvals(i)-1
        lambound = LVS(find(CNL(:)==j,1));
        vertsE(2*j-1:2*j,:) = [lambound 0; lambound 2];
        facesE(j-1,:) = [2*j-3 2*j-2 2*j 2*j-1];
        colsE(2*j-1:2*j) = [1; 1]*(j-2)/(Nvals(i)-3);
        if j==floor(Nvals(i)/2)
            lamdiv = lambound;
            vertsE(end-1:end,:) = [4*lambound 0; 4*lambound 2];
            facesE(end,:) = [2*Nvals(i)-3 2*Nvals(i)-2 2*Nvals(i) 2*Nvals(i)-1];
            colsE(end-1:end) = [1; 1];
        end
    end
    vertsE(:,1) = vertsE(:,1)/lamdiv;
    vertsE(:,2) = i+(vertsE(:,2)-1)/2;
    
    verts = [verts; vertsE];
    faces = [faces; facesE+infill];
    col = [col; colsE];
    quants(i,:) = [vertsE(3,1) quantile(vertsE(:,1),3) vertsE(end-2,1)];
    lamdivvec(i) = lamdiv;
    infill = infill+2*Nvals(i);
end
hold on
patch('Faces', faces, 'Vertices', verts, 'FaceVertexCData', col, 'FaceColor', 'interp', 'EdgeColor', 'none','FaceAlpha',.5)
infill = 0;
for i = 1:length(lamhvals)
    plot([quants(i,1) quants(i,5)],[i i],'k.','MarkerSize',6)
    plot([quants(i,2) quants(i,2)],[i-0.5 i+0.5],'k--','LineWidth',0.5)
    plot([1 1],[i-0.5 i+0.5],'k-','LineWidth',0.5)
    plot([quants(i,4) quants(i,4)],[i-0.5 i+0.5],'k--','LineWidth',0.5)
    if sources(i) == 0
        marktype = 'o';
    elseif sources(i) == 1
        marktype = 's';
    end
    if mechs(i) == 1
        markcolor = 'b';
    elseif mechs(i) == 0
        markcolor = 'w';
    elseif mechs(i) == -1
        markcolor = 'r';
    end
    errorbar(lamhvals(i)/lamdivvec(i),i,errsvals(i)/lamdivvec(i),'horizontal',marktype,'Markersize',8,...
        'MarkerFaceColor',markcolor,'MarkerEdgeColor','k','Color','k','LineWidth',1,'CapSize',8)
    infill = infill+Nvals(i);
end
hold off
colormap(mymap)
xlabel('Profile length, $\hat{\lambda}/\hat{\lambda}_{50}$','Interpreter','latex')
xlim([0 4])
xticks([0 1 2 3 4])
ylim([0.5 length(lamhvals)+0.5])
set(gca,'Ytick',[])
set(gca,'Yticklabel',[])
box on
annotation('textbox',[.525 .37 .1 .1],'String','E','LineStyle','none',...
    'Interpreter','latex','FontSize',16)

set(f,'Units','Inches')
fpos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[fpos(3),fpos(4)])
print(gcf,'Fig2','-dpdf')
print(gcf,'Fig2','-dpng','-r900')