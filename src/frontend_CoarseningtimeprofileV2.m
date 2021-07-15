function frontend_CoarseningtimeprofileV2(varargin)
%% This code analyzes trajectories to get the time profile of number of 
%% domains and domain properties at a cutoff value of 40 muM
%% Input variables
%path to load .mat file from
fileloadpathname=varargin{1};
loadmatfilename=varargin{2};
%% Output variables
% Nclusters - number of clusters as a cell array. Nclusters{1} will contain
% the number of clusters found in snapshot at time=5s from all replicates.
% clustervol - volume of clusters detected in micron^3
% FracActin - fraction of actin found in the cluster
% Domainmass - mass of actin in cluster 
% Data organization of the last 3 variables are the same as Nclusters
% variable.
%Height of the reaction volume.
[~, b] = ismember('Zspan',varargin);
if(b) Zspan = str2num(varargin{b+1});
    disp(['Using Zspan ',varargin{b+1}]);captag='';
end
%Size of pixel in nm
pixelsize = [100 100 100];
%Volume of the pixel
pixelVinL = pixelsize(1)*pixelsize(2)*pixelsize(3)*1e-27*1e3;
pixelvolumeinmum3 = pixelsize(1)*pixelsize(2)*pixelsize(3)*1e-9;
disp(['Pixel size ',num2str(pixelsize)]);
Rxnvoldim = [2000 2000 Zspan];
%Number of pixels
Npixels = ceil(Rxnvoldim./pixelsize);
% Npixels(1)*Npixels(2)*Npixels(3)*64/1e9
% 100nm -> 204Kb
% 50nm -> 15.4 Mb
% 25nm -> 9.8 Mb
% 10nm -> 1.92 Gb
% 5nm -> 15.36 Gb
% 1nm -> 1920 Gb
%Initialize instances of class pixel
for i = 1:Npixels(1)*Npixels(2)*Npixels(3)
    p(i) = pixel;
end
%Set number of pixels
p.setgetNpixels(Npixels);
disp('pixel created');
% Intiialize property variables
cutoff = 40;% micromolar cutoff
Nclusters = cell(1,400);%Number of clusters 
clustervol = cell(1,400);%Volume of clusters in \mum^3
FracActin = cell(1,400);%Fraction of total actin in a cluster
Domainmass = cell(1,400);%Total actin in a cluster
%% Get all files
%Fetch filenames of all replicates
filelist = ls([fileloadpathname,'/',loadmatfilename,'_S*.mat']);
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
%%
resetneeded = true;
%% Go through filelist and compute clusters 
for fid = 1:numel(filelist)
    load(filelist{fid});   
    disp('matfile loaded');
    disp(filelist{fid});
    %Find last snapshot
    for snap = 1:size(r(1).s,2)
        filcoord = r(1).s(snap).f.coord_cell1;
        if(isempty(filcoord))
            break;
        end
    end
    timevec = floor(r(1).time_vector(1:snap));
    %Snapshots are analyzed every 5 seconds
    ttargetvec = 5:5:max(timevec);
    for ttargetidx = 1:numel(ttargetvec)
        ttarget = ttargetvec(ttargetidx);
        snaplist = find(timevec==ttarget);
        for iter = 1
            snap = snaplist(iter);
            filcoord = r(1).s(snap).f.coord_cell1;
            disp(['pixellating snap ',num2str(snaplist(iter))]);
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
            plotstatus = false;
            countoccupied = 0;
            for fpos = 1:size(filcoord,1)
                % interpolate all monomers in each cylinder of the filament
                fil = reshape(filcoord{fpos}',3,[])';
                Coord_cat = interpolateallmonomersinfilemant(fil);
                if(plotstatus)
                    tempcoord = Coord_cat./pixelsize;
                    plot3(tempcoord(:,3),tempcoord(:,2),tempcoord(:,1),'ro-','LineWidth',2);hold on;
                end
                binvec = ceil(Coord_cat./pixelsize);
                onedpos = (binvec(:,3)-1).*crosssection + (binvec(:,2)-1).*Npixels(1) + binvec(:,1);
                onedpos = onedpos';
                for x=onedpos
                    if(~p(x).occupiedstatus)
                        countoccupied = countoccupied + 1;
                    end
                    % deem the pixel to be occupied, accounted, and then
                    % count the number of F-actin monomers
                    p(x).occupiedstatus = true;
                    p(x).accountedstatus = false;
                    p(x).count = p(x).count + 1;
                end
            end
            disp('pixel occupancy updated');
            %% Each entry is 1 monomers long.
            VinL=pi*1e3*1e3*7.5e3*1e-27*1e3;
            monomerspercount = 1;
            factor = monomerspercount*1e6/(6.023e23*VinL);
            Nmontotal = sum([p(:).count]);
            Capprox = factor*Nmontotal;
            factorpixel = monomerspercount*1e6/(6.023e23*pixelVinL);
            %% Update mean concentration of F-actin in each pixel
            for i = 1:size(p,2)
                pvecids = [p.getneighbors(i),i];
                p(i).meanconc = mean([p(pvecids).count])*factorpixel;
            end
            %% Calculate properties
            for j = 1:numel(cutoff)
                disp(['Calculating for cutoff ',num2str(cutoff(j))]);
                pcopy = p;
                C=[pcopy(:).meanconc];
                unaccountedpixels =find([p(:).accountedstatus]==false);
                nonzeroelems = intersect(find(C>=cutoff(j)),unaccountedpixels);
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
                    for cid = 1:size(clusters,2)
                        currcluster = clusters{cid};
                        cvolvec(cid) = numel(currcluster)*pixelvolumeinmum3;
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
    save(['coarseningtimeprofileV2-',loadmatfilename], 'Nclusters', 'clustervol','FracActin','Domainmass');
end
end
%% AUXILLARY FUNCTIONS
%Function to traverse tree network
function clusters = traversetree(p, countoccupied)
clusters = {};
%get all the neighbors of all pixels in one cell array
countunaccountedcubes = countoccupied;
neighborcell = cell(numel(p),1);
for i = 1:numel(p)
    neighborcell{i} = pixel.getneighbors(i);
end
%Begin recursive search
while(countunaccountedcubes >0)
    %Choose an unaccounted pixel at random
    unaccountedcubes = find([p(:).accountedstatus]==false);
    xx=randi(numel(unaccountedcubes));
    startnode = unaccountedcubes(xx);
    %Recursively search through it's neighbors and keep searching till you
    %hit an accounted pixel
    clusterpidvec=[];
    [clusterpidvec, countunaccountedcubes, p] = ...
        recursivetreesearch(startnode, clusterpidvec, countunaccountedcubes, p, neighborcell);
    if(numel(clusterpidvec))
        clusters = [clusters, clusterpidvec];
    end
end
end
% Recursive search function
function [clusterpidvec, countunaccountedcubes, p] = recursivetreesearch(nodeid, ...
    clusterpidvec, countunaccountedcubes,p,neighborcell)
%If it is unaccounted, consider it
if(p(nodeid).accountedstatus == false)
    p(nodeid).accountedstatus = true;
    countunaccountedcubes = countunaccountedcubes - 1;
    clusterpidvec = [clusterpidvec, nodeid];
    %Get the neighbors
    neighborlist = neighborcell{nodeid};
    npos  = [p(neighborlist).accountedstatus] == 0;
    neighborlist = neighborlist(npos);
    %Recursively search through neighbors till you hit an occupied pixel.
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