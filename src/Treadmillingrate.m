function [netminusendcell, netplusendcell, nfilcell, timecell] = Treadmillingrate(filepath,loadmatfilename)
% File calculates steady-state treadmilling rate based on MEDYAN trajectories.
%Load the monomer snapshot file
%% Input variables
% filepath - string - path to .mat files.
% loadmatfilename - string - name of the .mat file.
%% Output variables

% netminusendcell - cell array - stores the net minus end depolymerization
% between successive snapshots. netminusendcell{1} would be as long as
% sum_snap=N-x^N filaments(snap). So each entry coresponds to the change in
% number of monomers of a filament between successive snapshots.

% netplusendcell - cell array - stores the net plus end polymerizaton
% between successive snapshots. Data organization is same as
% netminusendcell

% nfilcell - cell array - stores the number of filaments. numel(nfilcell)
% gives the number of replicates, and nfilcell{1} gives the vector of
% filament counts in each snapshot considered.

%timecell - cell array - stores the time difference between successive
%snapshots. 
%%

%Fetch filenames of all replicates
filelist = ls([filepath,'/',loadmatfilename,'_MS*.mat']);
netminusendcell={};%Net loss in the minus end
netplusendcell={};% Net gain in the plus end
nfilcell={};% Number of filaments
timecell=cell(1,numel(filelist));%Time stamp of snapshots considered.
count = 1;
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
for fid = 1:numel(filelist)
    disp(filelist{fid});
    %load the trajectory .mat file
    load([filelist{fid}]);
    PPvec ={};%Plus end polymerizaton
    PMvec ={};%Minus end polymerization
    DPvec ={};%Plus end depolymerization
    DMvec ={};%Minus end depolymerization
    Factinvec ={};
    numFilvec = {};%Number of filaments
    % Starts from 1500th snapshot (in our trajectories for the current project,
    %we generate 1 snapshot per second)
    monomerfirstsnap = 1500;
    for runs = 1
        PPvecperrun =[];
        PMvecperrun =[];
        DPvecperrun =[];
        DMvecperrun =[];
        Factinvecperrun =[];
        numFilvecperrun = [];
        for msid = monomerfirstsnap:size(r(runs).ms,2)
            %If trajectory is long enough, get data
            if(msid>monomerfirstsnap)
                %Get the monomer snapshot of run and a snapshot
                monomers = r(runs).ms(msid).m;
                polyplus = sum(monomers.polyplusend);
                polyminus = sum(monomers.polyminusend);
                depolyplus = sum(monomers.depolyplusend);
                depolyminus = sum(monomers.depolyminusend);
                factin = sum(monomers.nummonomers);
                PPvecperrun = [PPvecperrun, polyplus];
                PMvecperrun = [PMvecperrun, polyminus];
                DPvecperrun = [DPvecperrun, depolyplus];
                DMvecperrun = [DMvecperrun, depolyminus];
                Factinvecperrun = [Factinvecperrun, factin];
                numFilvecperrun = [numFilvecperrun, numel(monomers.filID)];
            else
                disp('trajectory is too short');
            end
        end
        PPvec = [PPvec;PPvecperrun];
        PMvec = [PMvec;PMvecperrun];
        DPvec = [DPvec;DPvecperrun];
        DMvec = [DMvec;DMvecperrun];
        Factinvec = [Factinvec;Factinvecperrun];
        numFilvec = [numFilvec;numFilvecperrun];
        timecell{count} = r(1).mtime_vector;
        %Net depolymerization in minus end
        netminusendcell{count} = DMvec{1}-PMvec{1};
        %Net polymerizaton in plus end
        netplusendcell{count} = PPvec{1}-DPvec{1};
        nfilcell{count} = numFilvec;
        count = count + 1;
    end
end
end