load('JL_HVA_interneuronDensity_PV.mat')
load('JL_HVA_interneuronDensity_SOM.mat')

all_PV = logical(prod(cell2mat(cellfun(@(x) ~isnan(x(1,:)),density_PV,'un',0)')));
all_SOM = logical(prod(cell2mat(cellfun(@(x) ~isnan(x(1,:)),density_SOM,'un',0)')));
mean_PV_density = cell2mat(cellfun(@(x) nanmean(x(:,all_PV),2),density_PV,'un',0));
mean_SOM_density = cell2mat(cellfun(@(x) nanmean(x(:,all_SOM),2),density_SOM,'un',0));
area_names = {'LM','AL','PM','AM'};
layer_names = {'L2/3','L4','L5','L6'};
cmap = brewermap(4,'Set2');
cmap2 = flipud(brewermap(6,'Greys'));
%%

z_score_PV_density = reshape(zscore(mean_PV_density(:)),4,4);
z_score_SOM_density = reshape(zscore(mean_SOM_density(:)),4,4);

%%
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
%%

for area_i = 1:4
    temp_PV = z_score_PV_density(:,area_i);
    temp_SOM = z_score_SOM_density(:,area_i);
    interLayerDist(area_i,:) = pdist([temp_PV, temp_SOM]);
end     

figure;

fast_errbar(1:4,interLayerDist,2,'continuous',false,'stats',true)
xticks(1:4);xticklabels(area_names);xlim([0 5]);


%% compare interarea and interlayer distances
figure;
fast_errbar(1:2,[{interLayerDist(:)} {interAreaDist(:)}],1,'continuous',false,'stats',true);

%%
subplot(1,2,2); hold on;

for layer_i = 1:4
    plot([1 2 2],[interAreaDist(layer_i,1) interAreaDist(layer_i,2) interAreaDist(layer_i,3)],[marker_shape{layer_i}],'Color',(cmap(1,:)),'MarkerFaceColor',cmap(1,:));
    plot([1 2 2],[interAreaDist(layer_i,6) interAreaDist(layer_i,5) interAreaDist(layer_i,4)],[marker_shape{layer_i}],'Color',(cmap(4,:)),'MarkerFaceColor',cmap(4,:));
end
tempX = interAreaDist(:,[1 6]);
tempY = interAreaDist(:,[2 3 4 5]);

fast_errbar([1 2],[{tempX(:)} {tempY(:)}],1,'continuous',false);

xlim([0.5 2.5]); fix_axes(gcf,fs,'connection type','distance'); xticks([1 2]);
axis square
title(ttest2(tempX(:),tempY(:)))
