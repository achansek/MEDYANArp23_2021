function plottimeprofiledomains()
%% This version plots number of domains and concentration of actin in
%% domains as two different plots.
% Plots the number of domains and actin concentration in domains as time
% profiles. Please download shadedErrorBar function by Rob Campbell from
% https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
% 3. Please download colorSpectrum function by Ingo Hermann from
% https://www.mathworks.com/matlabcentral/fileexchange/57756-colorspectrum
colordelta = 1;
cutoff=27;
Nclustermean = []; Nclusterstd = [];
chordlenmean = []; chordlenstd = [];
Rgmean = []; Rgstd = [];
Fracmean = []; Fracstd = [];
FracDmassmean = []; FracDmassstd = [];
meancvol=[];stdcvol=[];
addpath E:\NINDS\April2017_unequaltrajectories\network_analysis\NINDS\Cyltry\Coarsening\timeprofile\Coarsening_matfiles_DT2;
errorbarvec=[];
TotalActin_mumoles=20*pi*1*1*7.5*1e-18*1e3;
color = colorSpectrum(6);
arptagcell = {'Arp1nM','Arp5nM','Arp10nM','Arp25nM','Arp50nM'};
enatagcell={'Ena1nM'};
extratag='Cyltry-';
kofftagcell={''};
volboxplot = [];
volboxids = [];
massboxplot = [];
massboxids = [];
nclusterboxplot = [];
nclusterids = [];
pdatavol1=[];pdatavol2=[];pdatavol3=[];pdatavol4=[];pdatavol5=[];pdatavol6=[];
pdatamass1=[];pdatamass2=[];pdatamass3=[];pdatamass4=[];pdatamass5=[];pdatamass6=[];
legendcell={'1','5','10','25','50'};
fig1 = figure('units','inch','position',[1,1,8,7],'Color','w');hax1=axes;
fig2 = figure('units','inch','position',[1,1,8,7],'Color','w');hax2=axes;
%%
Xaxisval = (5:5:2000)/1e3;
boxplotid = 1;
myosinalphatagcell = {'M-A-0-1-alpha-0-01'};
extratagcell={'Cyltry-'};
for matagid = 1:numel(myosinalphatagcell)
    matag = myosinalphatagcell{matagid};
    for koffid = 1:numel(kofftagcell)
        kofftag = [kofftagcell{koffid}];
        if(isempty(kofftagcell{koffid}))
            extratag='Cyltry-';
        else
            extratag='kinetics-';
            kofftag=[kofftag,'-'];
        end
        for etagid = 1:numel(extratagcell)
            extratag = extratagcell{etagid};
            for atagid = 1:numel(arptagcell)
                for enatagid = 1:numel(enatagcell)
                    arptag = [arptagcell{atagid}];
                    colorvec = color(atagid,:);
                    
                    if(~isempty(enatagcell{enatagid}))
                        arptag = [arptagcell{atagid},'-',enatagcell{enatagid}];
                    end
                    matloadfilepath = 'E:\NINDS\April2017_unequaltrajectories\network_analysis\NINDS\Cyltry\Coarsening\timeprofile\Coarsening_matfiles_DT2';
                    load([matloadfilepath,'\','coarseningtimeprofileV2-',extratag,matag,'-',kofftag,arptag,'.mat']);
                    Nfilteredclusters = cell(1,400);
                    filteredclustervols=cell(1,400);
                    lastsnap = 0;
                    for i = 1:400
                        X = clustervol{i};
                        n = Nclusters{i};
                        if(~isempty(n))
                            lastsnap = i;
                        end
                        nref = [0,cumsum(n)];
                        if(isempty(X))
                            disp(['Breaking at ',num2str(i)]);
                            Nfilteredclusters{i} = Nfilteredclusters{i-1};
                            continue;
                        end
                        ncorrected = [];
                        % Filter domains based on volume cutoff
                        locs=(find(X>cutoff*1e-3));
                        filteredclustervols{i} = X(locs);
                        for j = 2:numel(nref)
                            temp = numel(find(locs>nref(j-1)&locs<=nref(j)));
                            ncorrected = [ncorrected,temp];
                        end
                        Nfilteredclusters{i} = ncorrected;
                    end
                    %% moving average Number of domains
                    for i = 1:400
                        nclustersvec=[];
                        movdelta = 1;
                        for j = max(i-movdelta,1):min(400,i+movdelta)
                            nclustersvec = [nclustersvec, Nfilteredclusters{j}];
                        end
                        Nclustermean(i) = mean(nclustersvec);
                        Nclusterstd(i) = std(nclustersvec);
                    end
                    % last 100s
                    firstsnap = lastsnap - 20;
                    massvec=[];concvec=[];
                    clustervolvec=[];
                    clustermassvec=[];
                    nvec=[];
                    totalmass=294344;
                    for i=firstsnap:lastsnap
                        clustervolvec=[clustervolvec,clustervol{i}];
                        if(dontuse)
                            clustermassvec=[clustermassvec,Domainmass{i}];
                        end
                        nvec = [nvec, Nfilteredclusters{i}];
                    end
                    clustervolvec=clustervolvec./(100*100*100*1e-9);
                    locs = (find(clustervolvec>cutoff));
                    clustervolvec = clustervolvec(locs);
                    if(dontuse)
                        clustermassvec = clustermassvec(locs);
                    end
                    %Number of domains
                    axes(hax1);
                    shadedErrorBar(Xaxisval,Nclustermean,Nclusterstd,'lineProps',{'Color',colorvec,'LineWidth',2});hold on;
                    %% Concentration
                    massvec=[];concvec=[];
                    MdensityTP=zeros(1,400);
                    SdensityTP=zeros(1,400);
                    delta = 5;
                    for i=1:400
                        clustervolvec=[];
                        clustermassvec=[];
                        for j = max(1,i-delta):min(400,i+delta)
                            clustervolvec=[clustervolvec,clustervol{j}];
                            clustermassvec=[clustermassvec,Domainmass{j}];
                        end
                        clustervolvec=clustervolvec./(100*100*100*1e-9);
                        locs = (find(clustervolvec>cutoff));
                        clustervolvec = clustervolvec(locs);
                        if(dontuse)
                            clustermassvec = clustermassvec(locs);
                        end
                        density = 1e6*clustermassvec./(6.023e23*clustervolvec*1e6*1e-27*1e3);
                        MdensityTP(i) = mean(density);
                        SdensityTP(i) = std(density);
                    end
                    axes(hax2);
                    shadedErrorBar(Xaxisval,MdensityTP,SdensityTP,'lineProps',{'Color',colorvec,'LineWidth',2});hold on;
                end
            end
        end
    end
end
makepretty(2,legendcell);
ylim([40,350]);
ylabel({'Concentration of actin in','coarsened domains (\muM)'});
axes(hax1);
makepretty(1,legendcell);
ylim([0,12]);
ylabel({'Number of domains'});
%Save figures
savefig(fig1,'numdomains.fig');
print(fig1,'numdomains','-dpng','-r300')
savefig(fig2,'Acitnconcentration.fig');
print(fig2,'Actinconcentration','-dpng','-r300')
end
function makepretty(i,legendcell)
set(gca,'LineWidth',2);
set(gca,'FontSize',26);
xticklabels('auto');
yticklabels('auto');
leg = legend(legendcell,'NumColumns',1,'Location','NorthWest');
title(leg,'[Arp2/3] in nM');
xlabel('Time (\times 10^3 s)');
xlim([0,2]);
end
