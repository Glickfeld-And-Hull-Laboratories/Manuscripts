% load data 
load('JL_HVA_interneuronDensity_PV.mat')
load('JL_HVA_interneuronDensity_SOM.mat')

% plotting params
area_names = {'LM','AL','PM','AM'};
layer_names = {'L2/3','L4','L5','L6'};
cmap = brewermap(4,'Set2');
cmap2 = flipud(brewermap(6,'Greys'));
PV_color = [0.4 0.3 0.8];
SOM_color = [0.3 0.7 0.4];

% only take expt with all areas
all_PV = logical(prod(cell2mat(cellfun(@(x) ~isnan(x(1,:)),density_PV,'un',0)'))); 
all_SOM = logical(prod(cell2mat(cellfun(@(x) ~isnan(x(1,:)),density_SOM,'un',0)'))); 

% mean over animals 
mean_PV_density = cell2mat(cellfun(@(x) nanmean(x(:,all_PV),2),density_PV,'un',0));
mean_SOM_density = cell2mat(cellfun(@(x) nanmean(x(:,all_SOM),2),density_SOM,'un',0));

z_score_PV_density = reshape(zscore(mean_PV_density(:)),4,4);
z_score_SOM_density = reshape(zscore(mean_SOM_density(:)),4,4);
%% Figure 3B
figure;
ax(1) = subplot(1,4,1); hold on;
tempplot = cell2mat(cellfun(@(x) sum(x),density_PV,'un',0)');
tempplot = tempplot(:,all_PV)./max(tempplot(:,all_PV));
plotSpread(mat2cell(tempplot,ones(4,1),size(tempplot,2)),'distributionMarkers','o','distributionColor',PV_color);
fast_errbar(1:4,tempplot,2,'stats',true)
fix_axes(gcf,9,'Area','density');
axis square; 
ylim([0.5 1.1])

xticks(1:4);xticklabels(area_names);xlim([0.5 4.5]);

ax(2) = subplot(1,4,2); hold on;
tempplot = cell2mat(cellfun(@(x) sum(x),density_SOM,'un',0)');
tempplot = tempplot(:,all_SOM)./max(tempplot(:,all_SOM));
plotSpread(mat2cell(tempplot,ones(4,1),size(tempplot,2)),'distributionMarkers','o','distributionColor',SOM_color);
fast_errbar(1:4,tempplot,2,'stats',true)
fix_axes(gcf,9,'Area','density');
axis square; 
xticks(1:4);xticklabels(area_names);xlim([0.5 4.5]);


linkaxes(ax)
%% Figure 3C
% plot params
ylims = [-4 4];
fs = 10;

PV_nCellsByArea = repelem([{sum(all_PV)}],1,4);
SOM_nCellsByArea = repelem([{sum(all_SOM)}],1,4);

figure; 
for layer_i = 1:4
tmp_density_cell = cellfun(@(x) x(layer_i,all_PV)',density_PV,'un',0);
tmp_density_zscore = nanscore(cell2mat(tmp_density_cell'));
tmp_density_cell = mat2cell(tmp_density_zscore,cell2mat(PV_nCellsByArea));

subplot(2,4,layer_i); hold on;
plotSpread(tmp_density_cell,'spreadWidth',1.5,'DistributionColors',PV_color,'DistributionMarkers','o','binWidth',1);
fast_errbar(1:4,tmp_density_cell',1,'Continuous',false);

title(num2str(layer_names{layer_i}));
fix_axes(gcf,fs,'Area','z-scored cell density'); axis square;
ylim(ylims);
xticks(1:4);xticklabels(area_names);xlim([0.5 4.5]);

tmp_density_cell = cellfun(@(x) x(layer_i,all_SOM)',density_SOM,'un',0);
tmp_density_zscore = nanscore(cell2mat(tmp_density_cell'));
tmp_density_cell = mat2cell(tmp_density_zscore,cell2mat(SOM_nCellsByArea));

subplot(2,4,layer_i+4); hold on;
plotSpread(tmp_density_cell,'spreadWidth',1.5,'DistributionColors',SOM_color,'DistributionMarkers','o','binWidth',1);
fast_errbar(1:4,tmp_density_cell',1,'Continuous',false);

title(num2str(layer_names{layer_i}));
fix_axes(gcf,fs,'Area','z-scored cell density'); axis square;
ylim(ylims);
xticks(1:4);xticklabels(area_names);xlim([0.5 4.5]);

end
%% Figure 3D-E
fs = 20;
figure;

marker_shape = {'o','s','^','d'};
subplot(1,2,1); hold on;
for area_i = 1:4
    for layer_i = 1:4
    plot(z_score_PV_density(layer_i,area_i),z_score_SOM_density(layer_i,area_i),marker_shape{layer_i},'MarkerSize',10,'Color',cmap(area_i,:),'LineWidth',1.1,'MarkerFaceColor',cmap(area_i,:))
    end
end
axis square; matchxy('min',false);
fix_axes(gcf,fs,'Norm. PV Density','Norm. SOM Density');
legend(layer_names);

subplot(1,2,2); hold on

for layer_i = 1:4

    temp_PV = z_score_PV_density(layer_i,:)';
    temp_SOM = z_score_SOM_density(layer_i,:)';
    interAreaDist(layer_i,:) = pdist([temp_PV, temp_SOM]);
end

subplot(1,2,2);

for layer_i = 1:4
    plot([1 2 3],[interAreaDist(layer_i,1) interAreaDist(layer_i,2) interAreaDist(layer_i,3)],[marker_shape{layer_i}], 'MarkerSize',15,'MarkerFaceColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5]);
    plot([1 2 3],[interAreaDist(layer_i,6) interAreaDist(layer_i,5) interAreaDist(layer_i,4)],[marker_shape{layer_i}],'MarkerSize',15,'MarkerFaceColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5]);
end
tempX = interAreaDist(:,[1 6]);
tempY = interAreaDist(:,[2 5]);
tempZ = interAreaDist(:,[3 4]);

fast_errbar([1 2 3],[{tempX(:)} {tempY(:)} {tempZ(:)}],1,'continuous',false,'SD',true,'stats',true);

ylim([0 1.5]);xlim([0.5 3.5]); fix_axes(gcf,fs,'connection type','distance'); xticks([1 2 3]);
xticklabels([{'A-P'},{'M-L'},{'A-P + M-L'}]');
title([ttest2(tempX(:),tempY(:)),ttest2(tempX(:),tempZ(:)),ttest2(tempY(:),tempZ(:))])


%% InterLayer Dist by Area
for area_i = 1:4
    temp_PV = z_score_PV_density(:,area_i);
    temp_SOM = z_score_SOM_density(:,area_i);
    interLayerDist(area_i,:) = pdist([temp_PV, temp_SOM]);
end     
figure;
fast_errbar(1:4,interLayerDist,2,'continuous',false)
fix_axes(gcf,10,'Area','Z dist');
xticks(1:4);xticklabels(area_names);xlim([0 5]);


%% InterLayer vs InterArea distance + stats
figure;
fast_errbar(1:2,[{interLayerDist(:)} {interAreaDist(:)}],1,'continuous',false,'stats',true);
fix_axes(gcf,10,'Distance type','Z dist');
xticks(1:2);xticklabels({'interLayer','interArea'});xlim([0 3]);

%% Supplemental Figure 1

ylims = [0 15];
figure; 

for layer_i = 1:4
    
tmp_density_cell = {};
tmp_density_cell = cellfun(@(x) x(layer_i,logical(all_PV))',density_PV,'un',0);
tmp_density_cell = cellfun(@(x) (x/1000),tmp_density_cell,'un',0)'; % just divide for easier plot lims
subplot(2,4,layer_i); hold on;
plot(cell2mat(tmp_density_cell')','Color',PV_color)
plot(mean(cell2mat(tmp_density_cell')),'Color','k')
fix_axes(gcf,fs,'Area','Raw cell density (x 10^3)'); axis square;
xticks(1:4);xticklabels(area_names);xlim([0.5 4.5]);
ylim(ylims);
title(num2str(layer_names{layer_i}));

tmp_density_cell = {};
tmp_density_cell = cellfun(@(x) x(layer_i,logical(all_SOM))',density_SOM,'un',0);
tmp_density_cell = cellfun(@(x) (x/1000),tmp_density_cell,'un',0)'; % just divide for easier plot lims
subplot(2,4,layer_i+4); hold on;
plot(cell2mat(tmp_density_cell')','Color',SOM_color)
plot(mean(cell2mat(tmp_density_cell')),'Color','k')
fix_axes(gcf,fs,'Area','Raw cell density (x 10^3)'); axis square;
ylim(ylims);
xticks(1:4);xticklabels(area_names);xlim([0.5 4.5]);

end
