clear all
close all
clc

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
    [~,insort] = sort(lamdata(in));
    lamsort(ins+1:ins+length(in)) = in(insort);
    ins = ins+length(in);
end
lamhvals = lamdata(lamsort);
errsvals = errdata(lamsort);
Nvals = Ndata(lamsort);
sources = source(lamsort);
mechs = mech(lamsort);

phivals = zeros(size(lamhvals));
Cs = zeros(1,sum(Ndata));
mlabels = cell(sum(Ndata),1);
xmax = 5*ceil(max(lamhvals)/5);
verts = zeros(2*(sum(Nvals)+2*length(Nvals)),2);
faces = zeros(sum(Nvals)+length(Nvals),4);
col = zeros(2*(sum(Nvals)+2*length(Nvals)),1);
infill = 0;
for i = 1:length(lamhvals)
    phivals(i) = fminsearch(@(phi) (DTCC(Nvals(i),lamhvals(i),phi)-0.95)^2,-0.1);
    Cstemp = zeros(1,Nvals(i));
    for j = 1:Nvals(i)
        Cstemp(j) = fminsearch(@(lam) (Cvals(Nvals(i),lam,phivals(i),j)-1)^2,lamhvals(i));
        Cs(infill+j) = Cstemp(j);
        mlabels{infill+j} = mlabeldata{lamsort(i)};
    end
    quants = quantile(Cstemp,3);
    stin = 2*(infill+2*(i-1));
    verts(stin+1:stin+2*(Nvals(i)+2),:) = [[0 Cstemp/quants(2) 2.5 0 Cstemp/quants(2) 2.5]' ...
        [(i-0.5)*ones(Nvals(i)+2,1); (i+0.5)*ones(Nvals(i)+2,1)]];
    col(stin+1:stin+2*(Nvals(i)+2),:) = [0 0:1/(Nvals(i)-1):1 1 0 0:1/(Nvals(i)-1):1 1]';
    faces(infill+i:infill+i+Nvals(i),:) = [(stin+1:stin+Nvals(i)+1)' (stin+2:stin+Nvals(i)+2)' ...
        (stin+Nvals(i)+4:stin+2*Nvals(i)+4)' (stin+Nvals(i)+3:stin+2*Nvals(i)+3)'];    
    infill = infill+Nvals(i);
end

figure(1)
figmult = 0.5;
sperc = 1/3;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

s1 = subplot(221);
s1pos = get(s1,'Position');
set(s1,'Position',[s1pos(1) s1pos(2)+sperc*s1pos(4) s1pos(3) (1-sperc)*s1pos(4)])
lamhex = 10:1:40;
SP = 1:3;
rhojvals = zeros(length(SP),length(lamhex));
for i = 1:length(lamhex)
    disp(i/length(lamhex))
    for j = SP
        rhojvals(j,i) = genrhoj(50,100,lamhex(i),-0.05,j);
    end
end
semilogy(lamhex,rhojvals)
hold on
plot([lamhex(1) lamhex(end)],[1 1],'k--')
hold off
legend({'1D','2D','3D'},'Location','NorthEast','Interpreter','latex')
legend boxoff
xlabel('Profile Length, $\hat{\lambda}=\lambda/a$','Interpreter','latex')
ylabel('$\rho_{j}$','Interpreter','latex')

s2 = subplot(222);
s2pos = get(s2,'Position');
set(s2,'Position',[s2pos(1) s2pos(2)+sperc*s2pos(4) s2pos(3) (1-sperc)*s2pos(4)])
rhoc = zeros(1,100);
lamopt = lamhex(find((rhojvals(1,:)-1).^2==min((rhojvals(1,:)-1).^2),1));
for i = 1:length(rhoc)
    rhoc(i) = fminsearch(@(lam) (Cvals(100,lam,-0.05,i)-1)^2,lamopt);
end
quants = quantile(rhoc,3);
rhoc = rhoc/quants(2);
quants = quants/quants(2);
minlamc = 0;
maxlamc = 2;
verts2 = zeros(2*length(rhoc)+4,2);
verts2(1:2,:) = [minlamc 0; minlamc 2];
verts2(end-1:end,:) = [maxlamc 0; maxlamc 2];
faces2 = zeros(length(rhoc)+1,4);
faces2(1,:) = [1 2 4 3];
col2 = zeros(2*length(rhoc)+4,1);
col2(1:2) = [0; 0];
col2(end-1:end) = [1; 1];
for i = 1:length(rhoc)
    verts2(2*i+1:2*i+2,:) = [rhoc(i) 0; rhoc(i) 2];
    faces2(i+1,:) = [2*i+1 2*i+2 2*i+4 2*i+3];
    col2(2*i+1:2*i+2) = [1; 1]*(i-1)/(length(rhoc)-1);
end
hold on
patch('Faces', faces2, 'Vertices', verts2, 'FaceVertexCData', col2, 'FaceColor', 'interp', 'EdgeColor', 'none','FaceAlpha',.5)
plot([rhoc(1) rhoc(end)],[1 1],'k.','MarkerSize',6*figmult)
plot([quants(1) quants(1)],[0 2],'k--','LineWidth',0.5*figmult)
plot([quants(2) quants(2)],[0 2],'k-','LineWidth',0.5*figmult)
plot([quants(3) quants(3)],[0 2],'k--','LineWidth',0.5*figmult)
hold off
colormap(mymap)
ax1 = gca;
ax1.XLim = [0 2];
ax1.XTick = [0 1 2];
xlabel('Profile Length, $\hat{\lambda}/\hat{\lambda}_{50}$','Interpreter','latex')
ax1.YTick = [];
ax2 = axes('Position',ax1.Position,'Color','none');
ax2.XAxisLocation = 'top';
ax2.XLim = ax1.XLim;
ax2.XTick = [rhoc(1),quants(1),quants(2),quants(3),rhoc(end)];
ax2.XTickLabel = {'0%','25%','50%','75%','100%'};
ax2.YTick = [];

s3 = subplot(223);
s3pos = get(s3,'Position');
set(s3,'Position',[s3pos(1) s3pos(2) s3pos(3) (1+sperc)*s3pos(4)])
hold on
p1 = plot(NaN,NaN,'o','Markersize',9*figmult,'MarkerFaceColor','k','MarkerEdgeColor','k');
p2 = plot(NaN,NaN,'s','Markersize',9*figmult,'MarkerFaceColor','k','MarkerEdgeColor','k');
p3 = plot(NaN,NaN,'o','Markersize',9*figmult,'MarkerFaceColor','r','MarkerEdgeColor','k');
p4 = plot(NaN,NaN,'o','Markersize',9*figmult,'MarkerFaceColor','b','MarkerEdgeColor','k');
p5 = plot(NaN,NaN,'o','Markersize',9*figmult,'MarkerFaceColor','w','MarkerEdgeColor','k');
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
    errorbar(lamvals(lamsort(i)),i,errvals(lamsort(i)),'horizontal',marktype,'Markersize',12*figmult,...
        'MarkerFaceColor',markcolor,'MarkerEdgeColor','k','Color','k','LineWidth',figmult,'CapSize',12*figmult)
end
hold off
legend([p1 p2 p3 p4 p5],{'{\it Drosophila}','Zebrafish','Evidence for DT','Evidence for SDC',...
    sprintf('Evidence for\nmultiple mechanisms')},'Interpreter','latex','Location','SouthEast','FontSize',6)
legend boxoff
xlabel('Profile length, $\lambda$ ($\mu$m)','Interpreter','latex')
ylim([0.5 length(lamhvals)+0.5])
yticks(1:length(lamhvals))
yticklabels(mlabeldata(lamsort))
set(gca,'TickLabelInterpreter','latex')

s4 = subplot(224);
s4pos = get(s4,'Position');
set(s4,'Position',[s4pos(1) s4pos(2) s4pos(3) (1+sperc)*s4pos(4)])
hold on
patch('Faces', faces, 'Vertices', verts, 'FaceVertexCData', col, 'FaceColor', 'interp', 'EdgeColor', 'none','FaceAlpha',.5)
infill = 0;
for i = 1:length(lamhvals)
    quants = quantile(Cs(infill+1:infill+Nvals(i)),3);
    plot([Cs(infill+1) Cs(infill+Nvals(i))]/quants(2),[i i],'k.','MarkerSize',6*figmult)
    plot([quants(1) quants(1)]/quants(2),[i-0.5 i+0.5],'k--','LineWidth',0.5*figmult)
    plot([quants(2) quants(2)]/quants(2),[i-0.5 i+0.5],'k-','LineWidth',0.5*figmult)
    plot([quants(3) quants(3)]/quants(2),[i-0.5 i+0.5],'k--','LineWidth',0.5*figmult)
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
    errorbar(lamhvals(i)/quants(2),i,errsvals(i)/quants(2),'horizontal',marktype,'Markersize',12*figmult,...
        'MarkerFaceColor',markcolor,'MarkerEdgeColor','k','Color','k','LineWidth',figmult,'CapSize',12*figmult)
    infill = infill+Nvals(i);
end
hold off
colormap(mymap)
xlabel('Profile length, $\hat{\lambda}/\hat{\lambda}_{50}$','Interpreter','latex')
xlim([0 2.5])
xticks([0 1 2])
ylim([0.5 length(lamhvals)+0.5])
set(gca,'Ytick',[])
set(gca,'Yticklabel',[])

print(gcf,'Fig2','-dpdf')
print(gcf,'Fig2','-dpng','-r900')