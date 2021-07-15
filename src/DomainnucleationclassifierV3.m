function [Eventcountercell,dNclustercell]=...
    DomainnucleationclassifierV3(startsnap,endsnap)
%Nsplit Nmerge Nnucleated Ndestructed quantified in a more strict fashion
%% Input variables 
% startsnap, endsnap, the first and last snaps to consider when generating
% the event tally
%% Output variables
%Eventcountercell holds the tally of different events in the time dframe
%dNclustercell - number of cluster details
%deltaNcell - tallies the different events that lead to change in number of
%domains between frames.
%Nparentcell - number of parents
%Fetch filenames of all coarsening time profile files
filelist = ls('coarseningtimeprofile*.mat');
%If Windows
filelisttemp={};
if ispc
    for fid = 1:size(filelist,1)
        if(numel(filelist(fid,:))==0)
            %if filename is empty, skip.
            continue;
        end
        filelisttemp = [filelisttemp, [filepath,filelist(fid,:)]];
    end
    filelist = filelisttemp;
end
kmat=[];
rsqmat=[];
Nsample=2000;
Eventcountercell={};
dNclustercell={};
Nparentcell={};
%each row represents the cluster number. first element is the number of
%times the count decreased by 1 and the second element is the number of
%times it stayed the same and the third element is the number of times
%the domain count decreased by 1.
for fid = 1:numel(filelist)
    Eventcounter=zeros(15,3);
    dNcluster=zeros(15,4);
    parentmasscell=cell(15,4);
    parentvolumecell=cell(15,4);
    offspringmasscell=cell(15,4);
    offspringvolumecell=cell(15,4);
    allparentmasscell=cell(15,1);
    allparentvolumecell=cell(15,1);
    % Load file
    load(filenamecell{fid});
    %Set threshold values
    nucleationthreshold = 0.99;
    volumethreshold = 27;
    if(size(contribcellofcells,2)>Nsample)
        temp=[contribcellofcells(1:Nsample)];
        for i = 1:size(contribcellofcells,2)/Nsample - 1
            temp=[temp;contribcellofcells(i*Nsample+1:(i+1)*Nsample)];
        end
        contribcellofcells=temp;
    end
    dNoverall=[0 0 0 0];
    Nparents=[];
    firstsnapcluster=0;
    deltaNcell={};
    count = 1;
    startrun=1;
    if(fid==5)
        startrun=1;
    end
    runcounter = 1;
    deltaNeventseries=[];
    for runid = startrun+2
        contribcell=contribcellofcells(runid,:);
        contribmatprev = contribcell{startsnap-1};
        parentsizevec = sum(contribmatprev,2);
        %Find domains that are large enough to be considered
        parentlocs = find(parentsizevec>volumethreshold);
        %Get volume
        parentvolume = parentsizevec(parentlocs);
        %Get actin mass
        parentmass = Domainmass{startsnap-1}(parentlocs);
        predicted=[];
        actual=[];
        %Go through time
        for snap = startsnap:endsnap
            deltaN=[0 0 0 0];
            %Contribmat stores the information about current snap domains
            %with reference to domains in the previous volume. Thus, we
            %can track the pixels in space as they get assigned and
            %reassigned between different domains as a function of time.
            contribmat = contribcell{snap};
            contribmatcpy = contribmat;
            if(isempty(contribmat)||size(contribmat,2)-1~=size(contribmatprev,1))
                continue;
            end
            %We refer to domains from the previous frame as parents and the
            %domains from the next frame as offsprings.
            parentcontribvec = sum(contribmat(:,1:end-1),1);
            offspringsizevec = sum(contribmat,2);
            % parse through contribmat
            offspringlocs = find(offspringsizevec>volumethreshold);
            offspringvolume=offspringsizevec(offspringlocs);
            offspringmass=Domainmass{snap}(offspringlocs);
            contribmat = contribmat(:,[parentlocs;size(contribmat,2)]);
            contribmat = contribmat(offspringlocs,:);
            %Find the previous frame domains that are not assigned to any domain
            %in this frame.
            unaccounted_parents = 1:size(contribmat,2)-1;
            unaccounted_offsprings = 1:size(contribmat,1);
            delta = numel(unaccounted_parents)-numel(unaccounted_offsprings);
            nnparents = numel(unaccounted_parents);
            if(delta<0)
                Eventcounter(nnparents,3) = Eventcounter(nnparents,3) + abs(delta);
            elseif(delta==0)
                Eventcounter(nnparents,2) = Eventcounter(nnparents,2) + 1;
            else
                Eventcounter(nnparents,1) = Eventcounter(nnparents,1) + abs(delta);
            end
            
            %Ndestructed - if parent contributes to none of the current domains
            locs = find(sum(contribmat(:,1:end-1),1)==0);
            unaccounted_parents = setdiff(unaccounted_parents,locs);
            deltaN(4) = deltaN(4)+numel(locs);
            if(numel(locs))
                dNcluster(numel(parentlocs),4)=dNcluster(numel(parentlocs),4)+numel(locs);
                parentmasscell{numel(parentlocs),4}=[parentmasscell{numel(parentlocs),4},parentmass(locs)];
                parentvolumecell{numel(parentlocs),4}=[parentvolumecell{numel(parentlocs),4},parentvolume(locs)'];
            end
            %Nnucleated - receives most volume from bulk
            nuclocs=[];
            for i = 1:size(contribmat,1)
                temp = contribmat(i,:);
                %Nnucleated - bulk alone
                fracbulk = temp(end)./sum(temp);
                if(fracbulk>=nucleationthreshold)
                    deltaN(3) = deltaN(3) + 1;
                    dNcluster(numel(parentlocs),3)=dNcluster(numel(parentlocs),3)+1;
                    nuclocs=[nuclocs,i];
                end
            end
            if(~isempty(nuclocs))
                offspringmasscell{numel(parentlocs),3}=...
                    [offspringmasscell{numel(parentlocs),3},offspringmass(nuclocs)];
                offspringvolumecell{numel(parentlocs),3}=...
                    [offspringvolumecell{numel(parentlocs),3},offspringvolume(nuclocs)'];
            end
            unaccounted_offsprings = setdiff(unaccounted_offsprings,nuclocs);
            delta = numel(unaccounted_offsprings)-numel(unaccounted_parents);
            contribmatSM = contribmat(unaccounted_offsprings,:);
            contribmatSM = contribmatSM(:,unaccounted_parents);
            if(delta<0)
                %Merge
                deltaN(2) = deltaN(2) + 1;
                dNcluster(numel(parentlocs),2)=dNcluster(numel(parentlocs),2)+1;
                parentmasscell{numel(parentlocs),2}=...
                    [parentmasscell{numel(parentlocs),2},parentmass(unaccounted_parents)];
                parentvolumecell{numel(parentlocs),2}=...
                    [parentvolumecell{numel(parentlocs),2},parentvolume(unaccounted_parents)'];
                offspringmasscell{numel(parentlocs),2}=...
                    [offspringmasscell{numel(parentlocs),2},offspringmass(unaccounted_offsprings)];
                offspringvolumecell{numel(parentlocs),2}=...
                    [offspringvolumecell{numel(parentlocs),2},offspringvolume(unaccounted_offsprings)'];
            elseif(delta>0)
                %Split
                deltaN(1) = deltaN(1) + 1;
                dNcluster(numel(parentlocs),1)=dNcluster(numel(parentlocs),1)+1;
                parentmasscell{numel(parentlocs),1}=...
                    [parentmasscell{numel(parentlocs),1},parentmass(unaccounted_parents)];
                parentvolumecell{numel(parentlocs),1}=...
                    [parentvolumecell{numel(parentlocs),1},parentvolume(unaccounted_parents)'];
                offspringmasscell{numel(parentlocs),1}=...
                    [offspringmasscell{numel(parentlocs),1},offspringmass(unaccounted_offsprings)];
                offspringvolumecell{numel(parentlocs),1}=...
                    [offspringvolumecell{numel(parentlocs),1},offspringvolume(unaccounted_offsprings)'];
            end
            Nparents(count) = numel(parentlocs);
            if(snap==startsnap)
                firstsnapcluster = firstsnapcluster+ numel(parentlocs);
            end
            dNoverall = dNoverall + deltaN;
            deltaNeventseries=[deltaNeventseries;deltaN];
            count = count + 1;
            allparentmasscell{numel(parentlocs),1}=...
                [allparentmasscell{numel(parentlocs),1},parentmass];
            allparentvolumecell{numel(parentlocs),1}=...
                [allparentvolumecell{numel(parentlocs),1},parentvolume'];
            contribmatprev = contribmatcpy;
            parentsizevec = offspringsizevec;
            parentlocs = offspringlocs;
            parentvolume = offspringvolume;
            parentmass = offspringmass;
        end
        runcounter=runcounter+1;
    end
    deltaNcell=[deltaNcell,deltaNeventseries];
    Nparentcell=[Nparentcell,Nparents];
    Eventcountercell=[Eventcountercell,Eventcounter];
    dNclustercell=[dNclustercell,dNcluster];
    disp(num2str([firstsnapcluster,dNoverall]));
end
end