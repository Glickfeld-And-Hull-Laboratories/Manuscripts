% For recreating Figure 7
% Requires JL2021_loadPair_connectivity.m and
% JL2021_pairStim_connectivityParams.m

%% analyze days and save data grouped
your_read_path = 'Z:\home\jen\Notebook files\HVA recording\'; % where the xlsx notebook file is 
your_save_path = 'Z:\home\jen\Analysis\HVA recording\Pair Connectivity\'; % where to save data

% JL2021_loadPair_connectivity(your_read_path,your_save_path);
%% load data and set plotting params
ccc;
your_read_path = 'Z:\home\jen\Notebook files\HVA recording\'; % where the xlsx notebook file is 
your_save_path = 'Z:\home\jen\Analysis\HVA recording\Pair Connectivity\'; % where to save data

PV = load([your_save_path,'PV_grouped.mat']);
SOM = load([your_save_path,'SOM_grouped.mat']);

cmap = brewermap(5,'Set2');
PV_color = [0.3 0.4 0.7];
SOM_color = [0 0.5 0.3];
area_names = {'LM';'AL';'PM';'AM'};
cell_names = {'PV';'SOM'};

%% Figure 7B-C plot amplitudes of connected by area
doStats = true;
% PV 
figure; subplot(1,2,1);hold on;
cellfun(@(x,y) plot(x,abs(y),'o','Color',PV_color),{1;2;3;4},PV.IN_amp_group','un',0)
fast_errbar(1:4,cellfun(@abs,PV.IN_amp_group,'un',0),1,'stats',doStats);
fix_axes(gcf,15,'Area','Pyr->PV first pulse amplitude');
axis square
xticklabels(area_names);xlim([0 5]);ylim([0 150]);

subplot(1,2,2); hold on;
cellfun(@(x,y) plot(x,y,'o','Color',PV_color),{1;2;3;4},PV.pyr_amp_group','un',0)
fast_errbar(1:4,PV.pyr_amp_group,1,'stats',doStats);
fix_axes(gcf,15,'Area','PV->Pyr first pulse amplitude');
axis square
xticklabels(area_names);xlim([0 5]);ylim([0 150]);

% SOM
figure;
subplot(1,2,1); hold on;
cellfun(@(x,y) plot(x,abs(y),'o','Color',SOM_color),{1;2;3;4},SOM.pyr_amp_group','un',0)
fast_errbar(1:4,cellfun(@abs,SOM.IN_amp_group,'un',0),1,'stats',doStats);
fix_axes(gcf,15,'Area','Pyr->SOM first pulse amplitude');
axis square
xticklabels(area_names);xlim([0 5]);ylim([0 150]);

subplot(1,2,2); hold on;
cellfun(@(x,y) plot(x,y,'o','Color',SOM_color),{1;2;3;4},SOM.pyr_amp_group','un',0)
fast_errbar(1:4,SOM.pyr_amp_group,1,'stats',doStats);
fix_axes(gcf,15,'Area','SOM->Pyr first pulse amplitude');
axis square
xticklabels(area_names);xlim([0 5]);ylim([0 150]);


%% figure 7D colormaps 
SOM_cmap = brewermap(200,'Greens');
PV_cmap = brewermap(200,'Blues');

figure;
ax(1) = subplot(1,2,1);
groupData = [SOM.INpyr_byArea; SOM.pyrIN_byArea; SOM.recip_byArea];
imagesc(groupData,[0 1]);
colormap(ax(1),SOM_cmap);
xticks(1:4)
xticklabels(area_names)
yticks(1:3)
yticklabels([{'SOM->Pyr'};{'Pyr->SOM'};{'SOM<->Pyr'}])
[XX,YY] = meshgrid(1:size(groupData,2),1:size(groupData,1));
text(XX(:),YY(:),num2str(round(groupData(:),2)),'HorizontalAlignment','center','FontSize',15);
fix_axes(gcf,10,'Area','Connection type'); axis square;
colorbar


ax(2) = subplot(1,2,2);
groupData = [PV.INpyr_byArea; PV.pyrIN_byArea; PV.recip_byArea];
imagesc(groupData,[0 1])
colormap(ax(2),PV_cmap);
xticks(1:4)
xticklabels(area_names)
yticks(1:3)
yticklabels([{'PV->Pyr'};{'Pyr->PV'};{'PV<->Pyr'}])
[XX,YY] = meshgrid(1:size(groupData,2),1:size(groupData,1));
text(XX(:),YY(:),num2str(round(groupData(:),2)),'HorizontalAlignment','center','FontSize',15);
fix_axes(gcf,10,'Area','Connection type'); axis square;
colorbar

%% Table 2 P1/P2 and P1/P10 statistics

doStats = true;
fs = 10;
figure;
subplot(2,2,1)
fast_errbar(1:2,[{cell2mat(PV.IN_P1P2_group')} {cell2mat(SOM.IN_P1P2_group')}],1,'continuous',false,'stats',doStats)
title('Pyr->IN P2/P1');
xlim([0.5 2.5]); xticks([1:2]);xticklabels(cell_names);
fix_axes(gcf,fs,'Cell type','Ratio')
subplot(2,2,2)
fast_errbar(1:2,[{cell2mat(PV.pyr_P1P2_group')} {cell2mat(SOM.pyr_P1P2_group')}],1,'continuous',false,'stats',doStats)
title('IN->Pyr P2/P1');
xlim([0.5 2.5]); xticks([1:2]);xticklabels(cell_names);
fix_axes(gcf,fs,'Cell type','Ratio')
subplot(2,2,3)
fast_errbar(1:2,[{cell2mat(PV.IN_P1P10_group')} {cell2mat(SOM.IN_P1P10_group')}],1,'continuous',false,'stats',doStats)
title('Pyr->IN P10/P1');
xlim([0.5 2.5]); xticks([1:2]);xticklabels(cell_names);
fix_axes(gcf,fs,'Cell type','Ratio')
subplot(2,2,4)
fast_errbar(1:2,[{cell2mat(PV.pyr_P1P10_group')} {cell2mat(SOM.pyr_P1P10_group')}],1,'continuous',false,'stats',doStats)
title('IN->Pyr P10/P1');
xlim([0.5 2.5]); xticks([1:2]);xticklabels(cell_names);
fix_axes(gcf,fs,'Cell type','Ratio')

%% Table 2 cell type first pulse amplitudes
doStats = true;
group_amp_IN = [{cell2mat(PV.IN_amp_group')} {cell2mat(SOM.IN_amp_group')}];
group_amp_pyr = [{cell2mat(PV.pyr_amp_group')} {cell2mat(SOM.pyr_amp_group')}];
figure;
subplot(1,2,1);
fast_errbar(1:2,group_amp_IN,1,'continuous',false,'stats',doStats)
xlim([0.5 2.5]); xticks([1:2]);xticklabels(cell_names);
fix_axes(gcf,fs,'Cell type','pA')
subplot(1,2,2);
fast_errbar(1:2,group_amp_pyr,1,'continuous',false,'stats',doStats)
xlim([0.5 2.5]); xticks([1:2]);xticklabels(cell_names);
fix_axes(gcf,fs,'Cell type','pA')