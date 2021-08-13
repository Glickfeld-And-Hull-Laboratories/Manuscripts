%% Supplementary figure 1B-D

ccc;
root_path = 'Z:\home\jen\Imaging\EMX axon imaging\To process\';
save_path = 'Z:\home\jen\Imaging\EMX axon imaging\Processed\';
redoAnalysis = false;
if redoAnalysis
    doSingleAxonAnalysis(root_path,save_path);
end
load([save_path,'single_axon_F.mat']);
animal_names = {'EMXCre402';'EMXCre403';'EMXCreM'};
area_names = {'LM';'AL';'PM';'AM'};

cmap = flipud(brewermap(5,'Set3'));
cmap_animal = flipud(brewermap(5,'Greys'));

figure;
subplot(1,2,1); hold on;
for animal_i = 1
    plot(group_hist.(animal_names{animal_i}).histogram/max(group_hist.(animal_names{animal_i}).histogram),'Color',[0.9 0.9 0.9]);
    plot(count_fit{animal_i}(1:255)/max(count_fit{animal_i}(1:255)),'Color',cmap(animal_i,:));
    plot([thresholdF_grouped(animal_i) thresholdF_grouped(animal_i)],[0 1],'--','Color',cmap(animal_i,:));
end
fix_axes(gcf,20,'pixel intensity','probability');
axis square

subplot(1,2,2);hold on;
for animal_i = 1:3
    plot(1:4,pixelF_grouped(animal_i,:),'-','Color',cmap_animal(animal_i,:));
    plot(0:5,repmat(thresholdF_grouped(animal_i),1,6),'--','Color',cmap_animal(animal_i,:));
end
fast_errbar(1:4,pixelF_grouped,1,'stats',true);
xticks(1:4); xticklabels(area_names);xlim([0.5 4.5]);
fix_axes(gcf,20,'area','avg. pixel intensity');
ylim([0 100]); 
axis square;



%%

function doSingleAxonAnalysis(root_path,save_path)
animal_names = {'EMXCre402';'EMXCre403';'EMXCreM'};
area_names = {'LM';'AL';'PM';'AM'};

addpath(genpath(root_path));
for animal_i = 1:numel(animal_names)
    file_names = vertcat(dir([root_path,animal_names{animal_i},'\*.tif']).name);
    nFiles = size(file_names,1);
    group_hist.(animal_names{animal_i}).histogram = zeros(1,255);
    for file_i = 1:nFiles
        load_file = file_names(file_i,:);
        file_ids = strsplit(load_file,'_');
        area_id = file_ids{3};
        
        tiff_raw = reshape(cell2mat(arrayfun(@(x) imread(load_file,x),1:31,'un',0)),[1024 1024 31 3]);
        tiff_max = max(tiff_raw(:,:,:,1),[],3);
        tiff_max(tiff_max == 0) = NaN;
        
        if isfield(group_hist.(animal_names{animal_i}),([area_id,'_maxImg']))
            group_hist.(animal_names{animal_i}).([area_id,'_maxImg']) = [group_hist.(animal_names{animal_i}).([area_id,'_maxImg']) {tiff_max}];
        else
            group_hist.(animal_names{animal_i}).([area_id,'_maxImg']) = {tiff_max};
        end
        count = histcounts(tiff_max,1:256);
        group_hist.(animal_names{animal_i}).histogram = sum([group_hist.(animal_names{animal_i}).histogram; count]);
    end
    count_fit{animal_i} = fit((1:255)',group_hist.(animal_names{animal_i}).histogram','Gauss1');
    [max_count,max_hist] = max(count_fit{animal_i}(1:255));
    fluo_width = max_hist-find(count_fit{animal_i}(1:255)>0.5*max_count,1,'first');
    threshold_val.(animal_names{animal_i}) = max_hist+3*fluo_width;
end


%
for animal_i = 1:numel(animal_names)
    for area_i = 1:4
        thresh_img = group_hist.(animal_names{animal_i}).([area_names{area_i},'_maxImg']);
        nROIs = numel(thresh_img);
        for roi_i = 1:nROIs
            thresh_pixels = thresh_img{roi_i}>threshold_val.(animal_names{animal_i});
            thresh_img{roi_i}(~thresh_pixels) = 0;
            group_hist.(animal_names{animal_i}).([area_names{area_i},'_maxImg']){roi_i} = thresh_img{roi_i};
            avgPixelF.(animal_names{animal_i}).([area_names{area_i},'_avgF'])(roi_i) = sum(thresh_img{roi_i}(:))/sum(thresh_pixels(:));
            avgPixelF.(animal_names{animal_i}).([area_names{area_i},'_nPixelsF'])(roi_i) = sum(thresh_pixels(:));
        end
        avgPixelF.(animal_names{animal_i}).avgF_group(area_i) = mean(avgPixelF.(animal_names{animal_i}).([area_names{area_i},'_avgF']));
    end
end
thresholdF_grouped = cell2mat(struct2cell(threshold_val));
pixelF_grouped = cell2mat(struct2cell(structfun(@(x) x.avgF_group,avgPixelF,'un',0)));

save([save_path,'single_axon_F.mat'],'thresholdF_grouped','pixelF_grouped','group_hist','count_fit','threshold_val');
end