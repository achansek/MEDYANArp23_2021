function plotter_concsweepanalyses(varargin)
% This code plots the resulting analyses from frontend_concsweepanalyses
%% Input for the function
% varargin{1} - matfileloadpath - path to find .mat files in
% varargin{2} - filenamecell - cell consisting of filename strings that will help the user
% to load .mat files.
% For example: {'Arp1nM','Arp5nM','Arp10nM','Arp25nM','Arp50nM'};
%% Prerequisites
% 1. Plots require the use of tight_subplot. Please download from
% https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w
% 2. Please download shadedErrorBar function by Rob Campbell from
% https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
% 3. Please download colorSpectrum function by Ingo Hermann from
% https://www.mathworks.com/matlabcentral/fileexchange/57756-colorspectrum
matfileloadpath = varargin{1};
loadmatfilenamecell = varargin{2};
fig = figure('units','inch','position',[1,1,18,8],'Color','w');
% [yspacing xspacing], [bottommargin, topmargin], [leftmargin, rightmargin]
[ha, pos] = tight_subplot(1,2,[.05 .08],[.15 .02],[.07 .02]);
color = colorSpectrum(numel(loadmatfilenamecell));
additionaldatamat = [];
X = [1:1:100,110:10:300,350:50:700,800:100,1500];
count = 1;
legendcell={};
for fid = 1:numel(loadmatfilenamecell)
    filenamestring = loadmatfilenamecell{fid};
    load(['coarseningV2-',filenamestring,'.mat']);
    [Nclusters,clustervol,FracActin]=parseresult(Nclusters,clustervol,FracActin,Domainmass);
    %1
    Nmean=[];Nstd=[];
    for i =1:numel(Nclusters)
        Nmean(i) = mean(Nclusters{i});
        Nstd(i) = std(Nclusters{i});
    end
    additionaldatamat = [additionaldatamat;Nmean(3),Nmean(7)];
    %3
    NmeanFA=[];NstdFA=[];
    for i =1:numel(FracActin)
        f = FracActin{i};
        f = f(find(~isnan(f)));
        NmeanFA(i) = mean(f);
        NstdFA(i) = std(f);
    end
    legendcell=[legendcell,arptag];
    axes(ha(1));
    shadedErrorBar(X,Nmean,Nstd,'lineProps',{'Color',color(count,:),'LineWidth',2});hold on;
    axes(ha(2));
    shadedErrorBar(X,NmeanFA,NstdFA,'lineProps',{'Color',color(count,:),'LineWidth',2});hold on;
    count = count + 1;
end
axes(ha(1));
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
xticklabels('auto')
xticklabels('auto');
yticklabels('auto')
xlabel('Threshold for F-actin Concentration in voxel, C_{F} (\muM)');
ylabel({'Mean number of domains \langleN(C\geqC_{F})\rangle'});
ylim([0,10]);
plot([20,20],[0,50],'k--','LineWidth',2);
set(gca,'XScale', 'log');
leg = legend({'1','5','10','25','50'});
title(leg,'[Arp2/3] in nM');
%%
axes(ha(2));
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
xticklabels('auto')
xticklabels('auto');
yticklabels('auto')
xlabel('Threshold for F-actin Concentration in voxel, C_{F} (\muM)');
ylabel({'Fraction of Actin in domains \chi(C\geqC_{F})'});
ylim([0,1]);
% plot([0,500],[0.5 0.5],'k--','LineWidth',2);
% plot([51 51],[0,1],'k--','LineWidth',2);
% plot([141 141],[0,1],'k--','LineWidth',2);
set(gca,'XScale', 'log');
plot([20,20],[0,50],'k--','LineWidth',2);
%%
figure;
plot(additionaldatamat(:,1),'r','LineWidth',2);hold on;
plot(additionaldatamat(:,2),'b','LineWidth',2);
set(gca,'FontSize',14);
set(gca,'LineWidth',2);
xticklabels({'1','5','10','25','50'});
xlabel('[Arp2/3] at [Ena] = 50nM');
ylabel('Mean number of domains');
leg = legend({'6','31'});
title(leg,{'Cut-off','Concentration','\muM'});
%save figure
savefig(fig,'Concsweep_domainnumber_fracactin_100s.fig');
print(fig,'Concsweep_domainnumber_fracactin_100s','-dpng','-r300')
end
%% AUXILLARY FUNCTION
function [Nclustersp,clustervolp,FracActinp]=parseresult(Nclusters,clustervol,FracActin,Domainmass)
len = size(Nclusters{1},2);
threshold = 27;
TotalActin=294483;
pixelvolumeinmu3=0.1*0.1*0.1;
Nclustersp=cell(1,len);
clustervolp=cell(1,len);
FracActinp=cell(1,len);
for j=1:size(Nclusters,2)
    nvec=Nclusters{j};
    vvec=clustervol{j};
    fvec=FracActin{j};
    dmvec=Domainmass{j};
    if(isempty(nvec))
        break;
    end
    ncsum = [0,cumsum(nvec)];
    for iter=1:numel(ncsum)-1
        if(iter<97)
            delta = 3;
        else
            delta = 0;
        end
        for i = max(1,iter-delta):min(numel(ncsum),iter+delta)
            ids=ncsum(i)+1:ncsum(i+1);
            vsubset=vvec(ids);
            dmsubset=dmvec(ids);
            locs = find(vsubset>threshold*pixelvolumeinmu3);
            Nclustersp{iter}=[Nclustersp{iter},numel(locs)];
            clustervolp{iter}=[clustervolp{iter},vsubset(locs)];
            FracActinp{iter}=[FracActinp{iter},min(sum(dmsubset(locs))./TotalActin,1.0)];
        end
    end
end
end
