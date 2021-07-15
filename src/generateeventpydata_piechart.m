function generateeventpydata_piechart()
%% This function generates event tally for domain events such as nucleation,
%% destruction, merge, and split. 
%% Entire trajectory
[~,dNclustercell]=DomainnucleationclassifierV3(3,2000);
piechartmatrix=[];
for idx=1:5
    DNprobArp1=dNclustercell{idx};
    temp=sum(DNprobArp1);
    piechartmatrix=[piechartmatrix;temp(1),temp(3),temp(2),temp(4)];
end
save('pdata_piechart-3-2000.mat','piechartmatrix');

%% 3-500s
[~,dNclustercell]=DomainnucleationclassifierV3(3,2000);
piechartmatrix=[];
for idx=1:5
    DNprobArp1=dNclustercell{idx};
    temp=sum(DNprobArp1);
    piechartmatrix=[piechartmatrix;temp(1),temp(3),temp(2),temp(4)];
end
save('pdata_piechart-3-500.mat','piechartmatrix');

%% 501-1000s
[~,dNclustercell]=DomainnucleationclassifierV3(3,2000);
piechartmatrix=[];
for idx=1:5
    DNprobArp1=dNclustercell{idx};
    temp=sum(DNprobArp1);
    piechartmatrix=[piechartmatrix;temp(1),temp(3),temp(2),temp(4)];
end
save('pdata_piechart-501-1000.mat','piechartmatrix');

%% 1001-1500s
[~,dNclustercell]=DomainnucleationclassifierV3(3,2000);
piechartmatrix=[];
for idx=1:5
    DNprobArp1=dNclustercell{idx};
    temp=sum(DNprobArp1);
    piechartmatrix=[piechartmatrix;temp(1),temp(3),temp(2),temp(4)];
end
save('pdata_piechart-1001-1500.mat','piechartmatrix');

%% 1501-2000s
[~,dNclustercell]=DomainnucleationclassifierV3(3,2000);
piechartmatrix=[];
for idx=1:5
    DNprobArp1=dNclustercell{idx};
    temp=sum(DNprobArp1);
    piechartmatrix=[piechartmatrix;temp(1),temp(3),temp(2),temp(4)];
end
save('pdata_piechart-1501-2000.mat','piechartmatrix');
end