ccc;

% plotting params
all_hvas = {'lm','al','pm','am'};
all_layers = {'L1','L23','L4','L5','L6'};
cmap = brewermap(5,'Set2');
cmap([3,4],:) = cmap([4,3],:);

% load data
load('JL_HVA_axon_density.mat')

% reshape data to have matching depth axis (resample to fit mean)
for hva_i = 1:numel(all_hvas)
    for layer_i = 1:numel(all_layers)
        yy_pix(layer_i,hva_i) = nanmean(cell2mat(cellfun(@numel,hva_aggs.emx.(all_hvas{hva_i}).(all_layers{layer_i}),'un',0)));
    end
end
yy_pix = floor(nanmax(yy_pix,[],2));
hva_aggs.emx.grouped_norm_data = cell(1,4);
for hva_i = 1:numel(all_hvas)
    for layer_i = 1:numel(all_layers)
        temp_dat = hva_aggs.emx.(all_hvas{hva_i}).(all_layers{layer_i})(end-5:end);
        reshape_dat = cell2mat(cellfun(@(x) resample(x,yy_pix(layer_i),numel(x)),temp_dat,'un',0));
        hva_aggs.emx.grouped_norm_data{hva_i} = cat(2,hva_aggs.emx.grouped_norm_data{hva_i},reshape_dat);
    end
end
clear plot_y

% extract fluorescence mean and subset for only expt with all areas 
for hva_i = 1:numel(all_hvas)
if isfield(hva_aggs.emx,(all_hvas{hva_i}))
temp_density{hva_i} = squeeze(hva_aggs.emx.(all_hvas{hva_i}).f_per_pix);
else
temp_density{hva_i} = NaN;
end
end

hasData = cell2mat(cellfun(@(x) ~isnan(x),temp_density,'un',0)');    
allData = logical(prod(hasData,1));  
subsetData = cellfun(@(x) x(:,allData),temp_density,'un',0);


%% Figure 2B-C

temp_denom = max(cell2mat(subsetData'));
normAcrossArea = cellfun(@(x) x./temp_denom,subsetData,'un',0);
normAcrossMean = cell2mat(cellfun(@(x) nanmean(x,2),normAcrossArea,'un',0));

temp_mat = cell2mat(normAcrossArea');
area_label = cellfun(@(x,y) repmat(x,y,size(temp_mat,2)),[{1},{2},{3},{4}],[{5} {5} {5} {5}],'un',0);
area_label = vertcat(area_label{:});
layer_id = repmat([1:numel(all_layers)]',4,size(temp_mat,2));
[p_FI,~,stats_result_FI] = anovan(temp_mat(:),[{area_label(:)},{layer_id(:)}],'varnames',{'area','layer'},'model','interaction');


figure;
subplot(1,3,1);imagesc(normAcrossMean,[0 1]);
colormap(gray)
xticks(1:4);xticklabels(all_hvas); yticks(1:5);yticklabels(all_layers); axis square;
set(gca,'FontSize',15);
h = colorbar; h.Limits = [0 1];
xlabel('Area');ylabel('Depth');
clear plot_y
subplot(1,3,2)
fast_errbar(1:5,normAcrossArea,2,'Color',cmap,'cells_as_x',false);
fix_axes(gcf,15,'Area','norm. fluorescence');
xticks(1:5);xticklabels(all_layers);axis square; pbaspect([3 2 1]);

clear plot_y
subplot(1,3,3)
y_den = max(cell2mat(cellfun(@(x) max(mean(x,1)),hva_aggs.emx.grouped_norm_data,'un',0)));
for hva_i = 1:numel(all_hvas)
    for mouse_i = 1:size(hva_aggs.emx.grouped_norm_data{hva_i},1)
        plot_y(mouse_i,:) = smooth(hva_aggs.emx.grouped_norm_data{hva_i}(mouse_i,:),20);
    end
    fast_errbar([],plot_y./y_den,1,'shaded',true,'color',cmap(hva_i,:));
end
axis square;
fix_axes(gcf,15,'Depth','norm. fluorescence');xlim([0 max(cumsum(yy_pix))]);
pbaspect([3 2 1]);
xticks(cumsum(yy_pix)); xticklabels(all_layers); ylim([0 1.2]);yticks(0:0.5:1);

%% Figure 2D-E
normWithinArea = cellfun(@(x) x./max(x),temp_density,'un',0);
normWithinMean = cell2mat(cellfun(@(x) nanmean(x,2),normWithinArea,'un',0));

temp_mat = cell2mat(normWithinArea');
area_label = cellfun(@(x,y) repmat(x,y,size(temp_mat,2)),[{1},{2},{3},{4}],[{5} {5} {5} {5}],'un',0);
area_label = vertcat(area_label{:});
layer_id = repmat([1:numel(all_layers)]',4,size(temp_mat,2));
[p_FI,~,stats_result_FI] = anovan(temp_mat(:),[{area_label(:)},{layer_id(:)}],'varnames',{'area','layer'},'model','interaction');

figure;
subplot(1,3,1);imagesc(normWithinMean,[0 1]);
colormap(gray)
xticks(1:4);xticklabels(all_hvas); yticks(1:5);yticklabels(all_layers); axis square;
set(gca,'FontSize',15);
h = colorbar; h.Limits = [0 1];
xlabel('Area');ylabel('Depth');
clear plot_y
subplot(1,3,2)
fast_errbar(1:5,normWithinArea,2,'Color',cmap,'cells_as_x',false);
fix_axes(gcf,15,'Area','norm. fluorescence');
xticks(1:5);xticklabels(all_layers);axis square; pbaspect([3 2 1]);

clear plot_y
subplot(1,3,3)
y_den = cell2mat(cellfun(@(x) max(mean(x,1)),hva_aggs.emx.grouped_norm_data,'un',0));
for hva_i = 1:numel(all_hvas)
    for mouse_i = 1:size(hva_aggs.emx.grouped_norm_data{hva_i},1)
        plot_y(mouse_i,:) = smooth(hva_aggs.emx.grouped_norm_data{hva_i}(mouse_i,:),20);
    end
    fast_errbar([],plot_y./y_den(hva_i),1,'shaded',true,'color',cmap(hva_i,:));
end
axis square;
fix_axes(gcf,15,'Depth','norm. fluorescence');xlim([0 max(cumsum(yy_pix))]);
pbaspect([3 2 1]);
xticks(cumsum(yy_pix)); xticklabels(all_layers); ylim([0 1.2]);yticks(0:0.5:1);

%% stats for L2/3 vs all other layers; LM vs all other layers
 
temp_mat = cell2mat(normAcrossArea');
area_label = cellfun(@(x,y) repmat(x,y,size(temp_mat,2)),[{1},{2},{3},{4}],[{5} {5} {5} {5}],'un',0);
area_label = vertcat(area_label{:});
layer_id = repmat([1:numel(all_layers)]',4,size(temp_mat,2));
[p_FI,~,stats_result_FI] = anovan(temp_mat(:),[{area_label(:)},{layer_id(:)}],'varnames',{'area','layer'},'model','interaction','display','off');

%% Figure 2 - supplemental 1

 figure;
for layer_i = 1:numel(all_layers)
    subplot(1,numel(all_layers),layer_i)
    tempplot = cellfun(@(x) x(layer_i,:),subsetData,'un',0);
    plot(cell2mat(tempplot'),'Color',[0.7 0.7 0.7]);
    fast_errbar(1:4,tempplot,2,'continuous',false);
    xticks(1:4);xticklabels(all_hvas);
    axis square;
    xlim([0 5]);
    fix_axes(gcf,8,'Area','Density (avg fluo.)');
    title(all_layers{layer_i});
end
