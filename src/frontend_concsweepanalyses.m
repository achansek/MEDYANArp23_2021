function frontend_concsweepanalyses(varargin)
% This code computes the number of domains and fraction of actin in
% domains.
%Path to files with .mat files
fileloadpathname = varargin{1};
%cell of .mat filename identifier string 
% For example: {'Arp1nM'};
loadmatfilename = varargin{2};
Zspan=7500;%The height of reaction volume in nm.
pixelsize = [100 100 100];%pixel size in nm
pixelVinL = pixelsize(1)*pixelsize(2)*pixelsize(3)*1e-27*1e3;
pixelvolumeinmum3 = pixelsize(1)*pixelsize(2)*pixelsize(3)*1e-9;
disp(['Pixel size ',num2str(pixelsize)]);
Rxnvoldim = [2000 2000 Zspan];
%Total number of pixels to use
Npixels = ceil(Rxnvoldim./pixelsize);
% Npixels(1)*Npixels(2)*Npixels(3)*64/1e9
% 100nm -> 204Kb
% 50nm -> 15.4 Mb
% 25nm -> 9.8 Mb
% 10nm -> 1.92 Gb
% 5nm -> 15.36 Gb
% 1nm -> 1920 Gb
%Initiate members of class pixel in an array
for i = 1:Npixels(1)*Npixels(2)*Npixels(3)
    p(i) = pixel;
end
%Set total number of pixels
p.setgetNpixels(Npixels);
disp('pixel created');

% Intiialize property variables
%Cutoff concentration (in muM) to be considered
cutoff = [1:1:100,110:10:300,350:50:700,800:100,1500];
Nclusters = cell(1,numel(cutoff));
clustervol = cell(1,numel(cutoff));
FracActin = cell(1,numel(cutoff));
Domainmass = cell(1,numel(cutoff));
%%
resetneeded = true;
runid = 1;
%Fetch filenames of all replicates of the snapshot file
filelist = ls([filepath,'/',loadmatfilename,'_S*.mat']);
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
%Loop through all replicates.
for fid = 1:numel(filelist)
    disp(filelist{fid});
    %load the trajectory .mat file
    load([filelist{fid}]);
    disp('matfile loaded');
    %Find last snapshot
    for snap = 1:size(r(1).s,2)
        filcoord = r(1).s(snap).f.coord_cell1;
        if(isempty(filcoord))
            break;
        end
    end
    timevec = floor(r(1).time_vector(1:snap));
    %Last 100s from each replicate is considered
    ttargetvec = timevec(end)-100:10:timevec(end);
    for ttargetidx = 1:numel(ttargetvec)
        ttarget = ttargetvec(ttargetidx);
        snaplist = find(timevec==ttarget);
        for iter = 1
            snap = snaplist(iter);
            filcoord = r(1).s(snap).f.coord_cell1;
            disp(['pixellating run ',num2str(runid),' snap ',num2str(snaplist(iter))]);
            %% reset pixels
            if(resetneeded)
                for pidx = 1:numel(p)
                    p(pidx).occupiedstatus = false;
                    p(pidx).accountedstatus = false;
                    p(pidx).count = 0;
                    p(pidx).meanconc = 0;
                end
            end
            %%Mark cubes outside the boundary as occupied.
            p = pixelcylinderinitialize(p,pixelsize);
            disp('pixel initialized');
            crosssection = Npixels(1)*Npixels(2);
            % Go through filaments in this snap
            countoccupied = 0;
            for fpos = 1:size(filcoord,1)
                % interpolate cylinder to get coordinates of all monomers
                fil = reshape(filcoord{fpos}',3,[])';
                % check cube of each interpolated bead
                Coord_cat = interpolateallmonomersinfilemant(fil);
                binvec = ceil(Coord_cat./pixelsize);
                onedpos = (binvec(:,3)-1).*crosssection + (binvec(:,2)-1).*Npixels(1) + binvec(:,1);
                onedpos = onedpos';
                for x=onedpos
                    %If empty,mark the pixel as occupied
                    if(~p(x).occupiedstatus)
                        countoccupied = countoccupied + 1;
                    end
                    p(x).occupiedstatus = true;
                    p(x).accountedstatus = false;
                    p(x).count = p(x).count + 1;%Increase number of monomers by 1
                end
            end
            disp('pixel occupancy updated');
            %% Each entry is 1 monomers long.
            VinL=pi*1e3*1e3*7.5e3*1e-27*1e3;
            monomerspercount = 1;
            %factor to convert monomers to local concentration in muM
            factor = monomerspercount*1e6/(6.023e23*VinL);
            Nmontotal = sum([p(:).count]);
            Capprox = factor*Nmontotal;
            factorpixel = monomerspercount*1e6/(6.023e23*pixelVinL);
            %% Update mean concentration of F-actin in each pixel based on neighboring pixels
            for i = 1:size(p,2)
                pvecids = [p.getneighbors(i),i];
                p(i).meanconc = mean([p(pvecids).count])*factorpixel;
            end
            %% Calculate properties at each threshold concentration
            for j = 1:numel(cutoff)
                disp(['Calculating for cutoff ',num2str(cutoff(j))]);
                pcopy = p;
                C=[pcopy(:).meanconc];
                unaccountedpixels =find([p(:).accountedstatus]==false);
                nonzeroelems = intersect(find(C>cutoff(j)),unaccountedpixels);
                %Set rest of the pixels as accounted.
                zeroelems = setdiff(1:numel(p),nonzeroelems);
                for k = zeroelems
                    pcopy(k).accountedstatus = true;
                end
                countoccupied = numel(nonzeroelems);
                countemptycubes = Npixels(1)*Npixels(2)*Npixels(3) - countoccupied;
                disp(['Fraction unoccupied ',num2str(pixelsize(1)*pixelsize(2)*pixelsize(3)*countemptycubes./(Rxnvoldim(1)*Rxnvoldim(2)*Rxnvoldim(3)))]);
                % Go through recursively to get all necessary pore clusters
                clusters = traversetree(pcopy, countoccupied);
                % Run statistics on the pores found
                Nclusters{ttargetidx} = [Nclusters{ttargetidx}, numel(clusters)];
                cvolvec = zeros(1,size(clusters,2));
                cmassvec = zeros(1,size(clusters,2));
                if(isempty(clusters))
                    cvolvec = 0;
                    cmassvec = 0;
                else
                    %Loop through clusters found
                    for cid = 1:size(clusters,2)
                        currcluster = clusters{cid};
                        %Volume of the cluster in mum^3
                        cvolvec(cid) = numel(currcluster)*pixelvolumeinmum3;
                        %Mass of a cluster is the total actin in a pixel
                        cmassvec(cid) = sum([p(currcluster).count]);
                    end
                end
                clustervol{ttargetidx} = [clustervol{ttargetidx}, cvolvec];
                Domainmass{ttargetidx} = [Domainmass{ttargetidx}, cmassvec];
                aa = sum([pcopy(nonzeroelems).count])/sum([pcopy(:).count]);
                if(~isnan(aa))
                    FracActin{ttargetidx} = [FracActin{ttargetidx}, aa];
                end
            end
            resetneeded = true;
        end
    end
    save(['coarseningV2-',loadmatfilename], 'Nclusters', 'clustervol', 'clusterRg','FracActin','Domainmass');
end
end
%% ADDITIONAL FUNCTIONS
%Function to traverse tree network
function clusters = traversetree(p, countoccupied)
clusters = {};
%Count number of unaccounted pixels
countunaccountedcubes = countoccupied;
%get all the neighbors of all pixels in one cell array
neighborcell = cell(numel(p),1);
for i = 1:numel(p)
    neighborcell{i} = pixel.getneighbors(i);
end
while(countunaccountedcubes >0)
    unaccountedcubes = find([p(:).accountedstatus]==false);
    %Choose a random unaccounted cube
    xx=randi(numel(unaccountedcubes));
    startnode = unaccountedcubes(xx);
    clusterpidvec=[];
    %Use recursive search to determine clusters
    [clusterpidvec, countunaccountedcubes, p] = ...
        recursivetreesearch(startnode, clusterpidvec, countunaccountedcubes, p, neighborcell);
    if(numel(clusterpidvec))
        clusters = [clusters, clusterpidvec];
    end
end
end
%Recursive search algorithm
function [clusterpidvec, countunaccountedcubes, p] = recursivetreesearch(nodeid, ...
    clusterpidvec, countunaccountedcubes,p,neighborcell)
%Given a node id that is not accounted yet, we recursively call this
%function to go over the neighbors till we hit another accounted pixel
if(p(nodeid).accountedstatus == false)
    p(nodeid).accountedstatus = true;
    countunaccountedcubes = countunaccountedcubes - 1;
    clusterpidvec = [clusterpidvec, nodeid];
    neighborlist = neighborcell{nodeid};
    npos  = [p(neighborlist).accountedstatus] == 0;
    neighborlist = neighborlist(npos);
    for neighbornodeid = neighborlist
        [clusterpidvec, countunaccountedcubes, p] = recursivetreesearch(neighbornodeid, ...
            clusterpidvec, countunaccountedcubes,p,neighborcell);
    end
end
end
%Function plots the circular ends of a cylinder
function p = pixelcylinderinitialize(p,pixelsize)
cylindercenter = [1000 1000]./pixelsize(1:2);
Radius = 1000./pixelsize(1);
for pid = 1:size(p,2)
    offsetintegers = p.getidx3dfrom1d(pid);
    if(norm(offsetintegers(1:2)-cylindercenter-[0.5 0.5]) > Radius)
        p(pid).occupiedstatus = false;
        p(pid).accountedstatus = true;
    end
end
end