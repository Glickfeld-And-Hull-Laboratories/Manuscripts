% plot by layer, A+B/B, A-B/B, A/B

figure;
set(gcf, 'position', [223    31   573   754]);
areas = {'LM','AL','PM','AM'};
layerLabel = {'2/3', '4', '5', '6'};

raw_color = [0.5 0.7 0.9];
avg_color = [0.2 0.3 0.5];
for i_layer = 1:4
    tmp_density = [];
    tmp_volume = [];
    tmp_counts = [];
    for i_area = 1:numel(areas)
            tmp_density(i_area, :) = density{i_area}(i_layer,:); % [Nareas, Nmice]
            tmp_volume(i_area,:) = volume{i_area}(i_layer,:); % [Nareas, Nmice]
            tmp_counts(i_area,:) = counts{i_area}(i_layer,:); % [Nareas, Nmice]
    end
    
    ax1(i_layer) = subplot(4,3, (i_layer-1)*3 +1);
    hold on
    plot(tmp_density, 'linewidth', 0.15,'Color',raw_color)
    fast_errbar(1:4,tmp_density, 2,'Color',avg_color);
    xticks(1:4); xticklabels(areas); axis square;
    box off; axis tight;
    ylabel('density (raw)')
    xlabel('Brain Area')
    title(sprintf('Layer %s', layerLabel{i_layer}));
    
    ax2(i_layer) = subplot(4,3, (i_layer-1)*3 +2);
    hold on
    plot(tmp_density./tmp_density(1,:), 'linewidth', 0.15,'Color',raw_color)
    fast_errbar(1:4,tmp_density./tmp_density(1,:), 2,'Color',avg_color);
    box off; axis tight
    xticks(1:4); xticklabels(areas);
    ylabel('density (A/B)')
    xlabel('Brain Area')
    title(sprintf('Layer %s', layerLabel{i_layer}));
    
    ax3(i_layer) = subplot(4,3, (i_layer-1)*3 +3)
    hold on
    plot((tmp_density-tmp_density(1,:))./tmp_density(1,:), 'linewidth', 0.15,'Color',raw_color)
    fast_errbar(1:4,(tmp_density-tmp_density(1,:))./tmp_density(1,:), 2,'Color',avg_color);
    box off; axis tight
    xticks(1:4); xticklabels(areas);
    ylabel('density (A-B/B)')
    xlabel('Brain Area')
    title(sprintf('Layer %s', layerLabel{i_layer}));  
end
linkaxes(ax1);linkaxes(ax2);linkaxes(ax3);

%%
figure; 
clear tmp_density
cmap = brewermap(4,'Set2');
for i_area = 1:4
    hold on;
    fast_errbar(1:4,density_PV{i_area}, 2,'Color',cmap(i_area,:));
    box off; axis tight; axis square;
    xticks(1:4); xticklabels(areas);
    ylabel('density (raw)')
    xlabel('Brain Area')
    
end

legend(layerLabel);

%%
all_layers = {'L1','L2/3','L4','L5/6'};
cmap = brewermap(4,'Set2');
areas = {'LM','AL','PM','AM'};
figure; hold on;
for area_i = 1:numel(areas)
    plot(histcounts(popdat(areas{area_i}).cellDepth,linspace(0,1000,100)),'Color',cmap(area_i,:),'LineWidth',1.5);
end
xticks([0,10,35,60]);xticklabels(all_layers');
    

%%
all_layers = {'L2/3','L4','L5','L6'};
figure;hold on;
tmp_density = cellfun(@(x) x*.01,density_SOM,'un',0);
for area_i = 1:4
    fast_errbar(1:4,tmp_density{area_i},2,'continuous',true,'Color',cmap(area_i,:));
end
fix_axes(gcf,10,'Layer','Density (cells/volume unit *** need to look up)');
xticks(1:4); xticklabels(all_layers');

%% within layer 2/3 area comparison
cmap = brewermap(4,'Set2');
doSubset= true;
plotDots = false;
doZScore = false;
if doSubset
    all_areas = prod(cell2mat(cellfun(@(x) ~isnan(x(1,:)),density,'un',0)'));
    nCellsByArea = repelem([{sum(all_areas)}],1,4);
else
    all_areas = cellfun(@(x) ones(1,size(x,2)),density,'un',0);
    nCellsByArea = cellfun(@(x) size(x,2),density,'un',0);
end
area_names = {'LM','AL','PM','AM'};
fs = 8;
temp = cellfun(@(x,y) repmat(x,y,1),area_names,nCellsByArea,'un',0);
area_ids = vertcat(temp{:});
    
figure; 

for layer_i = 1:4
    
tmp_density_cell = {};
if doSubset
    tmp_density_cell = cellfun(@(x) x(layer_i,logical(all_areas))',density,'un',0);
else
    tmp_density_cell = cellfun(@(x) x(layer_i,:)',density,'un',0);
end
if doZScore
    tmp_density_zscore = nanscore(cell2mat(tmp_density_cell'));
    tmp_density_cell = mat2cell(tmp_density_zscore,[cell2mat(nCellsByArea)]);
    ylims = [-4 4];
else
    tmp_density_cell = cellfun(@(x) (x/1000),tmp_density_cell,'un',0)';
    ylims = [0 15];
end

if plotDots
    fs = 10;
    subplot(1,4,layer_i); hold on;
    plotSpread(tmp_density_cell,'spreadWidth',1.5,'DistributionColors',cmap,'DistributionMarkers','o','binWidth',1);
    fast_errbar(1:4,tmp_density_cell',1,'Continuous',false,'stats',true);
else
    subplot(1,5,layer_i); hold on;
    plot(cell2mat(tmp_density_cell')','Color',[0.7 0.7 0.7])
    plot(mean(cell2mat(tmp_density_cell')),'Color','k')
end

[pval,~,stats] = anova1(cell2mat(tmp_density_cell)',area_ids,'off');
% figure;multcompare(stats)
% title(num2str(pval));
title(num2str(layer_names{layer_i}));
% xticks(1:4);xticklabels(area_names);xlim([0.5 4.5]);
fix_axes(gcf,fs,'Area','z-scored cell density'); axis square;
ylim(ylims);
% for area_i = 1:4
%     subplot(4,5,(5*(layer_i-1))+1+area_i); hold on;
%     tmp_density = cellfun(@(x) x(layer_i,:)./density{area_i}(layer_i,:),density,'un',0);
%     for sub_area = 1:4
%         plot(sub_area,tmp_density{sub_area},'o','Color',cmap(sub_area,:));
%     end
%     
%     fast_errbar(1:4,tmp_density,2,'Continuous',false);
%     [pval,~] = anova1(cell2mat(tmp_density),area_ids,'off');
%     title(num2str(pval));
%     
% fix_axes(gcf,fs,'Area','Norm. cell density');  
xticks(1:4);xticklabels(area_names);xlim([0.5 4.5]);
% hline(1,'k--');
% ylim([0.2 2]); axis square;
% end
end


%%
