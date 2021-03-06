function frontend_gettreadmillingrate(varargin)
% This code calculates the steady-state treadmilling rate of the network from
% .mat files generated by readMEDYANtraj.
%% Input for the function
% varargin{1} - matfileloadpath - path to find .mat files in
% varargin{2} - loadmatfilenamecell - cell consisting of filename strings that will help the user
% to load .mat files.
% For example: {'Arp0nM','Arp500pM','Arp1nM','Arp5nM','Arp10nM','Arp25nM','Arp50nM'};
% varargin{3} - xticklabelcell - cell consisting of tags to use as labels
% on xaxis of the plot
matfileloadpath = varargin{1};
loadmatfilenamecell = varargin{2};
xticklabelcell=varargin{3};
Nelems=numel(loadmatfilenamecell);
TrateMinuscell=cell(1,Nelems);
TratePluscell=cell(1,Nelems);
Tratetimecell=cell(1,Nelems);
TrateNfilcell=cell(1,Nelems);
for idx = 1:Nelems
    loadmatfilename = loadmatfilenamecell{idx};
    disp(loadmatfilename);
    [netminusendcell, netplusendcell, nfilcell, timecell]=...
        Treadmillingrate(matfileloadpath, loadmatfilename);
    TrateMinuscell{idx} = netminusendcell;
    TratePluscell{idx} = netplusendcell;
    Tratetimecell{idx} = timecell;
    TrateNfilcell{idx} = nfilcell;
    save('TRatemat.mat','loadmatfilenamecell','TrateMinuscell',...
        'TratePluscell','Tratetimecell','TrateNfilcell');
end
plotfigureandsave(xticklabelcell);
end
%% Plotter
function plotfigureandsave(xticklabelcell)
load('TRatemat.mat');
%% plotter
fig = figure('units','inch','position',[1,1,8,8],'Color','w');
meanvec=[];
stdvec=[];
idmat=1:numel(loadmatfilenamecell);
counter = 1;
for idx = idmat
    netminusendcell = TrateMinuscell{idx};
    netplusendcell = TratePluscell{idx};
    nfilcell = TrateNfilcell{idx};
    Tratedata=[];
    for i= 1:numel(netminusendcell)
        nm = netminusendcell{i};
        if(numel(nm)<300)
            disp(['Error ',loadmatfilenamecell{idx},'_MS',num2str(i)]);
            continue
        end
        numFilvec = nfilcell{i}{1}(end-300:end);
        numFilvec=1;
        nm = netminusendcell{i}(end-300:end)./numFilvec;
        np = netplusendcell{i}(end-300:end)./numFilvec;
        Tratedata=[Tratedata,nm,np];
    end
    meanvec(counter) = mean(Tratedata);
    stdvec(counter) = std(Tratedata);
    counter = counter + 1;
end
bar(meanvec,'r');
hold on;
errorbar(meanvec,stdvec,'ko','LineWidth',2);
xticklabels(xticklabelcell);
xtickangle(90);
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
ylabel('Treadmilling rate (monomers/s)');
xlabel('[Arp2/3] in nM');
%save figure
savefig(fig,'TreadmillingRate.fig');
print(fig,'TreadmillingRate','-dpng','-r300')
end