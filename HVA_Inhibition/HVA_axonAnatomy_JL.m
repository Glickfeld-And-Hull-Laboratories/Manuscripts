ccc;
all_hvas = {'lm','al','pm','am'};
all_layers = {'L1','L23','L4','L5','L6'};
cmap = brewermap(5,'Set2');

load('JL_HVA_axon_density.mat')

%%

for hva_i = 1:numel(all_hvas)
    for layer_i = 1:numel(all_layers)
        yy_pix(layer_i,hva_i) = nanmean(cell2mat(cellfun(@numel,hva_aggs.emx.(all_hvas{hva_i}).(all_layers{layer_i}),'un',0)));
    end
end

yy_pix = floor(nanmax(yy_pix,[],2));

hva_aggs.emx.grouped_norm_data = cell(1,4);
for hva_i = 1:numel(all_hvas)
    for layer_i = 1:numel(all_layers)
        temp_dat = hva_aggs.emx.(all_hvas{hva_i}).(all_layers{layer_i});
        
        reshape_dat = cell2mat(cellfun(@(x) resample(x,yy_pix(layer_i),numel(x)),temp_dat,'un',0));
%         reshape_dat = reshape_dat(:,7:end-7);
        hva_aggs.emx.grouped_norm_data{hva_i} = cat(2,hva_aggs.emx.grouped_norm_data{hva_i},reshape_dat);
    end
end


clear plot_y


for hva_i = 1:numel(all_hvas)
if isfield(hva_aggs.emx,(all_hvas{hva_i}))
test{hva_i} = squeeze(hva_aggs.emx.(all_hvas{hva_i}).f_per_pix);
else
test{hva_i} = NaN;
end
end
%%
doZscore = false;
LMnorm = false;
L23norm = false;
maxNorm = true;
if LMnorm
    densityByArea = cellfun(@(x) x./test{1}(2,:),test,'un',0);
elseif L23norm
    densityByArea = cellfun(@(x) x./x(2,:),test,'un',0);
elseif maxNorm
    densityByArea = cellfun(@(x) x./max(x),test,'un',0);
else
    densityByArea = {};
    if doZscore
        for layer_i = 1:numel(all_layers)

            tmp_density_cell = {};
            tmp_density = cellfun(@(x) x(layer_i,:),test,'un',0);
            tmp_density_zscore = nanscore(cell2mat(tmp_density));
            temp_cell(layer_i,:,:) = reshape(tmp_density_zscore,11,4);
        end
    
        for area_i = 1:4
            densityByArea{area_i} = temp_cell(:,:,area_i);
        end
    else
        densityByArea = test;
    end
end

%%

temp_mat = cell2mat(densityByArea');
area_label = cellfun(@(x,y) repmat(x,y,size(temp_mat,2)),[{1},{2},{3},{4}],[{5} {5} {5} {5}],'un',0);
area_label = vertcat(area_label{:});
layer_id = repmat([1:numel(all_layers)]',4,size(temp_mat,2));
[p_FI,~,stats_result_FI] = anovan(temp_mat(:),[{area_label(:)},{layer_id(:)}],'varnames',{'area','layer'},'model','interaction');
%%
doSubset = true;

if doSubset
    hasData = cell2mat(cellfun(@(x) ~isnan(x),densityByArea,'un',0)');    
    allData = logical(prod(hasData,1));  
    subsetData = cellfun(@(x) x(:,allData),densityByArea,'un',0);
else
    subsetData = densityByArea;
end
figure;
for layer_i = 1:numel(all_layers)
    subplot(1,numel(all_layers),layer_i)
    tempplot = cellfun(@(x) x(layer_i,:),subsetData,'un',0);
%     plotSpread(tempplot,'distributionColors',[0.5 0.5 0.5],'distributionMarkers','o','binWidth',0.25);
    plot(cell2mat(tempplot'),'Color',[0.7 0.7 0.7]);
    fast_errbar(1:4,tempplot,2,'continuous',false,'stats',true);
    xticks(1:4);xticklabels(all_hvas);
    axis square;
%     ylim([0 1]);
    xlim([0 5]);
    fix_axes(gcf,8,'Area','Density (avg fluo.)');
    title(all_layers{layer_i});
end

%% norm to different layers
if L23norm
    figure;
    for layer_i = 1:numel(all_layers)
        densityByArea = cellfun(@(x) x./x(layer_i,:),test,'un',0);
        for layer_j = 1:numel(all_layers)
            subplot(numel(all_layers),numel(all_layers),layer_j+(layer_i-1)*numel(all_layers))
            fast_errbar(1:4,cellfun(@(x) x(layer_j,:),densityByArea,'un',0),2,'continuous',false,'stats',true);
            xticks(1:4);xticklabels(all_hvas);
            axis square;
    %     ylim([0 1]);
            xlim([0 5]);
            fix_axes(gcf,6,'Area',['Fluo (norm', all_layers{layer_i}]);
        end
    end
end

%%
test2 = cell2mat(cellfun(@(x) nanmean(x,2),subsetData,'un',0));
doNorm = false;
if doNorm
    temp_denom = max(test2(:));
    test2 = test2./temp_denom;
    test = cellfun(@(x) x./temp_denom,test,'un',0);
end
    
%%

figure;

subplot(1,3,1);imagesc(test2,[0 1]);
colormap(gray)
xticks(1:4);xticklabels(all_hvas); yticks(1:5);yticklabels(all_layers); axis square;
set(gca,'FontSize',15);
h = colorbar; h.Limits = [0 1];
xlabel('Area');ylabel('Depth');
clear plot_y
subplot(1,3,2)
fast_errbar(1:5,test,2,'Color',cmap);
fix_axes(gcf,15,'Area','norm. fluorescence');
xticks(1:5);xticklabels(all_layers);axis square; pbaspect([3 2 1]);yticks(1:4);
% ylim([0 0.004]);
ylim([0 1.2]);yticks(0:0.5:1.0);

clear plot_y
subplot(1,3,3)
if doNorm
    y_den = max(cell2mat(cellfun(@(x) max(mean(x,1)),hva_aggs.emx.grouped_norm_data,'un',0)));
else
    y_den = 1;
end
for hva_i = 1:numel(all_hvas)
    for mouse_i = 1:size(hva_aggs.emx.grouped_norm_data{hva_i},1)
        plot_y(mouse_i,:) = smooth(hva_aggs.emx.grouped_norm_data{hva_i}(mouse_i,:),20);
    end
    fast_errbar([],plot_y./y_den,1,'shaded',true,'color',cmap(hva_i,:));
end
axis square;
xticks(cumsum(yy_pix)); xticklabels(all_layers); ylim([0 1.2]);yticks(0:0.5:1);
fix_axes(gcf,15,'Depth','norm. fluorescence');xlim([0 max(cumsum(yy_pix))]);
pbaspect([3 2 1]);
% ylim([0 0.004]);
%%
