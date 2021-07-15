function runmeanalysisscript
%% Front end control panel file to analyze MEDYAN trajectories
fil_length = true;
treadmilling = true;
concsweepsnaps = true;
concsweepanalysis = true;
timeprofilesnaps = true;
timeprofileanalysis = true;
domainkineticeventclassifier = true;
%% Please edit the following two lines
inputcell{1} = 'path/to/file';
fnamecell = {'filename1','filename2','filename3','filename4'};

%% Compute filament length distribution and plot
if(fil_length)
    % Please ignore _MSrunid part of the filename when creating the cell
    % array.
    frontend_getfillength(inputcell);
end


%% Compute treadmilling rate and plot
if(treadmilling)
    %Edit the following line
    inputcell{3} = {'1', '2', '3', '4'};%Labels to use in X axis
    % Please ignore _MSrunid part of the filename when creating the cell
    % array.
    frontend_gettreadmillingrate(inputcell);
end


%% Compute the pixellated concentration from last snapshots and plot them at
% various cutoff concentrations
if(concsweepsnaps)
    frontend_getconcsweeplastsnap(inputcell);
end


%% Last 100s of trajectories were used to compute the properties of
% clusters at various cutoff concentration
%Note: uses results from timeprofileanalysis flag
if(concsweepanalysis)
    for idx = 1:numel(fnamecell)
        inputcell{2} = fnamecell;
        frontend_concsweepanalyses(inputcell);
    end
    inputcell{2} = fnamecell;
    plotter_concsweepanalyses(inputcell);
end


%% Time profile panel plots from changing Arp2/3
if(timeprofilesnaps)
    frontend_timeprofilesnaps(inputcell);
end


%% Time profile analysis of the clusters
if(timeprofileanalysis)
    for idx = 1:numel(fnamecell)
        inputcell{2} = fnamecell;
        frontend_CoarseningtimeprofileV2(inputcell)
    end
    frontend_clusterdata_last100s();
    
    
    %Run DomainVolume_violinplot.ipynb to get the violin plots
    py.DomainVolume_violinplot.ipynb
    %Run the following file to get the violon plots of domain number and
    %actin concentration in domains from the last 100s of trajectories at
    %all myosin concentrations studied
    py.Violinplots_myosin_domains_concentration.ipynb
    %Note: If this fails, run python code as a standalone code using your
    %choice of python IDE.
    
    %Plot time profile of number of domains and actin concentration
    plottimeprofiledomains();
end

%% Compute domain nucleation, destruction, domain split and merge events
if(domainkineticeventclassifier)
    generateeventpydata_piechart();
    
    %Run the folllowing file to get the pie charts of domain events.
    py.Domainevents_piecharts.ipynb
    %Note: If this fails, run python code as a standalone code using your
    %choice of python IDE.
end
end