function frontend_clusterdata_last100s()
%% This code computes the cluster properties such as cluster volume,
%% cluster count, actin in cluster, and actin concentration in cluster
%% from the last 100s of the trajectory
cutoff=27;% the clusters found must span atleast 27 pixels in volume
%20 micromolar actin converted to micromoles in the reaction volume or
%radius 1micron and height 7.5 micron
TotalActin_mumoles=20*pi*1*1*7.5*1e-18*1e3;
arptagcell={'Arp1nM','Arp5nM','Arp10nM','Arp25nM','Arp50nM'};
pdatanum1=[];pdatanum2=[];pdatanum3=[];pdatanum4=[];pdatanum5=[];pdatanum6=[];
pdatavol1=[];pdatavol2=[];pdatavol3=[];pdatavol4=[];pdatavol5=[];pdatavol6=[];
pdatamass1=[];pdatamass2=[];pdatamass3=[];pdatamass4=[];pdatamass5=[];pdatamass6=[];
pdataconc1=[];pdataconc2=[];pdataconc3=[];pdataconc4=[];pdataconc5=[];pdataconc6=[];
Ntotalcell=400;
Xaxisval = (5:5:2000)/60;
myosinalphatagcell = {'M-A-0-1-alpha-0-01','M-A-0-05-alpha-0-01','M-A-0-01-alpha-0-01'};
for matagid = 1:numel(myosinalphatagcell)
    matag = myosinalphatagcell{matagid};
    for atagid = 1:numel(arptagcell)
        arptag = arptagcell{atagid};
        load([matloadfilepath,'/','coarseningtimeprofile-',matag,'-',kofftag,arptag,'.mat']);
        Nfilteredclusters = cell(1,Ntotalcell);
        filteredclustervols=cell(1,Ntotalcell);
        lastsnap = 0;
        for i = 1:Ntotalcell
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
        %moving average
        for i = 1:Ntotalcell
            nclustersvec=[];
            movdelta = 0;%Right now moveing average is disabled as movdelta value is 0.
            for j = max(i-movdelta,1):min(Ntotalcell,i+movdelta)
                nclustersvec = [nclustersvec, Nfilteredclusters{j}];
            end
            Nclustermean(i) = mean(nclustersvec);
            Nclusterstd(i) = std(nclustersvec);
        end
        % last 100s, 5s at a time
        firstsnap = lastsnap - 20;
        massvec=[];concvec=[];
        clustervolvec=[];
        clustermassvec=[];
        clusternumvec=[];
        nvec=[];
        totalmass=294344;
        for i=firstsnap:lastsnap
            clustervolvec=[clustervolvec,clustervol{i}];
            clustermassvec=[clustermassvec,Domainmass{i}];
            nvec = [nvec, Nfilteredclusters{i}];
        end
        clustervolvec=clustervolvec./(100*100*100*1e-9);
        locs = (find(clustervolvec>cutoff));
        clustervolvec = clustervolvec(locs);
        clustermassvec = clustermassvec(locs);
        clusterconcvec =1e6*clustermassvec./(6.023e23*clustervolvec*1e6*1e-27*1e3);
        
        % python compatible mat-file generate
        if(atagid==1)
            pdatanum1 = [pdatanum1, nvec];
            pdatavol1 = [pdatavol1, clustervolvec./1e3];
            pdatamass1 = [pdatamass1, clustermassvec*100./totalmass];
            pdataconc1 = [pdataconc1, clusterconcvec];
        elseif(atagid==2)
            pdatanum2 = [pdatanum2, nvec];
            pdatavol2 = [pdatavol2, clustervolvec./1e3];
            pdatamass2 = [pdatamass2, clustermassvec*100./totalmass];
            pdataconc2 = [pdataconc2, clusterconcvec];
        elseif(atagid==3)
            pdatanum3 = [pdatanum3, nvec];
            pdatavol3 = [pdatavol3, clustervolvec./1e3];
            pdatamass3 = [pdatamass3, clustermassvec*100./totalmass];
            pdataconc3 = [pdataconc3, clusterconcvec];
        elseif(atagid==4)
            pdatanum4 = [pdatanum4, nvec];
            pdatavol4 = [pdatavol4, clustervolvec./1e3];
            pdatamass4 = [pdatamass4, clustermassvec*100./totalmass];
            pdataconc4 = [pdataconc4, clusterconcvec];
        elseif(atagid==5)
            pdatanum5 = [pdatanum5, nvec];
            pdatavol5 = [pdatavol5, clustervolvec./1e3];
            pdatamass5 = [pdatamass5, clustermassvec*100./totalmass];
            pdataconc5 = [pdataconc5, clusterconcvec];
        elseif(atagid==6)
            pdatanum6 = [pdatanum6, nvec];
            pdatavol6 = [pdatavol6, clustervolvec./1e3];
            pdatamass6 = [pdatamass6, clustermassvec*100./totalmass];
            pdataconc6 = [pdataconc6, clusterconcvec];
        end
    end
    if(matagid==1)
        save('pdataBoxplotV2-M-A-0-1.mat','pdata*');
    elseif(matagid==2)
        save('pdataBoxplotV2-M-A-0-05.mat','pdata*');
    else
        save('pdataBoxplotV2-M-A-0-01.mat','pdata*');
    end
end
disp('data generated! Run the python script now');
end