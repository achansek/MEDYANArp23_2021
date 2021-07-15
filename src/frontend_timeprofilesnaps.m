function frontend_timeprofilesnaps(varargin)
%% This code generates a panel plot at various time points along the trajectory
%% of actin density fields given a local actin concentration threshold of 40 micromolar.
fig = figure('units','inch','position',[1,1,11,1],'Color','w');
%nrows, ncols, [yspacing xspacing], [bottommargin, topmargin], [leftmargin, rightmargin]
[ha, pos] = tight_subplot(1,7,[.2 .05],[.07 .1],[.07 .07]);
pixelsize = [100 100 100];%pixel size in nm
%Volume of pixel
pixelVinL = pixelsize(1)*pixelsize(2)*pixelsize(3)*1e-27*1e3;
pixelvolumeinmum3 = pixelsize(1)*pixelsize(2)*pixelsize(3)*1e-9;
disp(['Pixel size ',num2str(pixelsize)]);
Rxnvoldim = [2000 2000 7500];%Reaction volume dimensions in nm
Npixels = ceil(Rxnvoldim./pixelsize);%Number of pixels in a frame
% Npixels(1)*Npixels(2)*Npixels(3)*64/1e9
% 100nm -> 204Kb
% 50nm -> 15.4 Mb
% 25nm -> 9.8 Mb
% 10nm -> 1.92 Gb
% 5nm -> 15.36 Gb
% 1nm -> 1920 Gb
%%
%Path to files with .mat files
fileloadpathname = varargin{1};
%cell of .mat filename identifier string
% For example: {'Arp1nM', 'Arp5nM', 'Arp10nM', 'Arp25nM', 'Arp50nM'};
loadmatfilenamecell = varargin{2};
for fid = 1:numel(loadmatfilenamecell)
    loadmatfilename = loadmatfilenamecell{fid};
    clear p;
    for i = 1:Npixels(1)*Npixels(2)*Npixels(3)
        p(i) = pixel;
    end
    p.setgetNpixels(Npixels);
    disp('pixel created');
    %%
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
    %Just use the first replicate run to generate the images
    for fid = 1
        disp(filelist{fid});
        %load the trajectory .mat file
        load([filelist{fid}]);
        disp('matfile loaded');
        runid = 1;
        
        % Get the filament coordinates
        for snap = 1:size(r(1).s,2)
            filcoord = r(1).s(snap).f.coord_cell1;
            if(isempty(filcoord))
                break;
            end
        end
        %time points along the time profile to be plotted
        snapref = [5 100 500 800 1200 1600 2000]+1;
        timevec = r(1).time_vector;
        %Find the corresponding snapshot
        snaplist = [];
        for i =1:numel(snapref)
            snaplist(i) = find(timevec<=snapref(i),1,'last');
        end
        count = 1;
        for snap = snaplist
            filcoordcell{count} = r(1).s(snap).f.coord_cell1;
            count = count + 1;
            resetneeded = true;
            for iter = 1:numel(snaplist)
                filcoord = filcoordcell{iter};
                disp(['pixellating run ',num2str(runid),' snap ',num2str(snaplist(iter))]);
                % reset pixels
                if(resetneeded)
                    for pidx = 1:numel(p)
                        p(pidx).occupiedstatus = false;
                        p(pidx).accountedstatus = false;
                        p(pidx).count = 0;
                        p(pidx).meanconc = 0;
                    end
                end
                % Mark cubes outside the boundary as occupied.
                p = pixelcylinderinitialize(p,pixelsize);
                disp('pixel initialized');
                crosssection = Npixels(1)*Npixels(2);
                % Go through filaments in this snap
                plotstatus = false;
                countoccupied = 0; 
                for fpos = 1:size(filcoord,1)
                    % interpolate all monomers in a filament
                    fil = reshape(filcoord{fpos}',3,[])';
                    Coord_cat = interpolateallmonomersinfilemant(fil);
                    if(plotstatus)
                        tempcoord = Coord_cat./1e3;
                        plot3(tempcoord(:,3),tempcoord(:,2),tempcoord(:,1),'ro-','LineWidth',2);hold on;
                    end
                    binvec = ceil(Coord_cat./pixelsize);
                    onedpos = (binvec(:,3)-1).*crosssection + (binvec(:,2)-1).*Npixels(1) + binvec(:,1);
                    onedpos = onedpos';
                    for x=onedpos
                        
                        if(~p(x).occupiedstatus)
                            countoccupied = countoccupied + 1;
                        end
                        p(x).occupiedstatus = true;
                        p(x).accountedstatus = false;
                        p(x).count = p(x).count + 1;
                    end
                end
                disp('pixel occupancy updated');
                % Each entry is 10 monomers long.
                VinL=pi*1e3*1e3*7.5e3*1e-27*1e3;
                monomerspercount = 1;
                factor = monomerspercount*1e6/(6.023e23*VinL);
                Nmontotal = sum([p(:).count]);
                Capprox = factor*Nmontotal;
                factorpixel = monomerspercount*1e6/(6.023e23*pixelVinL);
                % Update mean concentration of F-actin in each pixel
                for i = 1:size(p,2)
                    pvecids = [p.getneighbors(i),i];
                    p(i).meanconc = mean([p(pvecids).count])*factorpixel;
                    %         p(i).meanconc = mean(p(i).count)*factorpixel;
                end
                %
                cutoff = 40;%Generate time profiles at a threshold of 40 micromolar.
                C=[p(:).meanconc];
                Accounted = [p(:).accountedstatus];
                withincylinder = find(Accounted==false);
                for j = 1:numel(cutoff)
                    nonzeroelems = find(C>cutoff(j));
                    %remove out of boundary pixels
                    nonzeroelems = intersect(withincylinder,nonzeroelems);
                    Coords=zeros(numel(nonzeroelems),3);
                    count = 1;
                    for i = nonzeroelems
                        Coords(count,:) = pixelsize.*pixel.getidx3dfrom1d(i);
                        count = count + 1;
                    end
                    Coords = Coords./1e3;
                    plotid = fid*(numel(loadmatfilenamecell)-1)+j;
                    axes(ha(plotid));plot3(Coords(:,3),Coords(:,2),Coords(:,1),'rs','LineWidth',0.1,'MarkerFaceColor','r','MarkerSize',6);
                    tvec = r(1).time_vector;
                    title(num2str(round(tvec(snaplist(iter)),1)));
                    xlim([0,7.5]);
                    ylim([0,2]);
                    zlim([0,2]);
                    view([-10,12]);
                    set(gca,'LineWidth',2);
                    set(gca,'FontSize',20);
                    % Plot cylinder
                    plotcylinder();
                    set(gca,'visible','off')
                end
            end
        end
    end
end
%save figure
savefig(fig,'Timesweep_traj.fig');
print(fig,'Timesweep_traj','-dpng','-r300')
end
%% AUXILLARY FUNCTIONS
function p = pixelcylinderinitialize(p,pixelsize)
cylindercenter = [1000 1000]./pixelsize(1:2);
Radius = 1000./pixelsize(1);
count = 1;
for pid = 1:size(p,2)
    offsetintegers = p.getidx3dfrom1d(pid);
    if(norm(offsetintegers(1:2)-cylindercenter-[0.5 0.5]) > Radius)
        p(pid).occupiedstatus = false;
        p(pid).accountedstatus = true;
        count = count + 1;
    end
end
disp(['Number of pixels outside boundary=',num2str(count)]);
end
function plotcylinder()
hold on;
[X,Y,Z]=cylinder(1,50);
Z=Z*0.25;
surf(Z,Y+1,X+1,'FaceAlpha',0.4,'FaceColor','y','EdgeColor','k');
surf(Z+7.5-0.25,Y+1,X+1,'FaceAlpha',0.4,'FaceColor','y','EdgeColor','k');
hold off;
end