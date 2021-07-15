function  [Laxis, cdfmat, pdfmat, pdfmatCapped, pdfmatEna] = getfillengthCDFCap(matloadpath, matfilename)
%% Code clculates the filament length from number of monomers in each filament 
%% Input variables
% matloadpath - string - stores the path to the parsed .mat files.
% matfilename - string - stores the parsed monomer snapshot .mat file to
% load
%% Output variables
% Laxis - vector of bin mid points
% cdfmat - matrix of cumulative density function. Each row corresponds to a
% particular replicate.
% pdfmat - matrix of probability density function corresponding to filament 
% lengths. Each row corresponds to a particular replicate.
% pdfmatCapped - matrix of filament length probability density functions
% corresponding to Capped filament ends. Each row corresponds to a 
% particular replicate.
% pdfmatEna - matrix of filament length probability density functions
% corresponding to Ena bound filament ends. Each row corresponds to a 
% particular replicate.

%%
% Bins in increments of 50nm
Lbins = 0:50:1e4;
Laxis = 0.5*(Lbins(1:end-1) + Lbins(2:end));
cdfmat=[];
pdfmat = [];
pdfmatCapped=[];
pdfmatEna=[];
% Handles up to 10 replicates
for runid = 1:10
    % Load monomer snapshot file
    disp([matloadpath,'/',matfilename,'_MS',num2str(runid),'.mat']);
    if(isfile([matloadpath,'/',matfilename,'_MS',num2str(runid),'.mat']))
        rms = load([matloadpath,'/',matfilename,'_MS',num2str(runid),'.mat']);
        rms = rms.r;
        disp('loaded');
        %Takes the last 250s of snapshot
        minsnap = size(rms(1).ms,2)-250;
        %Skips if the trajectory is shorter than 250s.
        if(minsnap<0)
            disp('snap<250');
            continue;
        end
        L = zeros(1,numel(Laxis));
        LCapped = zeros(1,numel(Laxis));
        LEna = zeros(1,numel(Laxis));
        %Bins data from every 10 frames. Assuming 1s per frame and 250
        %final shaphots, this would mean that the snapshots that are 10s
        %apart are used.
        for snap = minsnap:10:size(rms(runid).ms,2)
            fillengthCap = [];
            fillengthEna = [];
            fillengthmatrix = 2.7*rms(runid).ms(snap).m.nummonomers;
            Plusendstatus = rms(runid).ms(snap).m.cap;
            %Bins the filament length
            L = L + histcounts(fillengthmatrix,Lbins);
            % If trajectory has internal inconsistencies, skip.
            % The number of plusends must match the number of filaments in
            % any given snapshot.
            if(numel(Plusendstatus)~=numel(fillengthmatrix))
                disp('skipping');
                continue
            else
                fillengthCap = fillengthmatrix(find(Plusendstatus==1));
                fillengthEna = fillengthmatrix(find(Plusendstatus==2));
                LCapped = LCapped + histcounts(fillengthCap ,Lbins);
                LEna = LEna + histcounts(fillengthEna ,Lbins);
            end
        end
        % Normalization factor
        normfactor = trapz(Laxis,L);
        % L is the probability density function
        L = L./normfactor;
        LCapped = LCapped./normfactor;
        LEna = LEna./normfactor;
        %calculate cumulative density function.
        cdfL = zeros(1,numel(Laxis));
        cdfL(1) = L(1);
        for k = 2:numel(Laxis)
            cdfL(k) = trapz(Laxis(1:k),L(1:k));
        end
        cdfmat = [cdfmat;cdfL];
        pdfmat = [pdfmat;L];
        pdfmatCapped = [pdfmatCapped;LCapped];
        pdfmatEna = [pdfmatEna;LEna];
        runid = runid + 1;
    else
        disp('not found');
    end
end
end
