% load data
ccc;
area_names = {'LM','AL','PM','AM'};
excel_path = 'grouped_interneuron_properties_ratio.xlsx';
summary_table = readtable(excel_path);

is_vClamp = cell2mat(cellfun(@(x) ~isempty(x),summary_table.HFS,'un',0));
filenames = cellfun(@(x,y) [x,y],summary_table.Date,summary_table.HFS,'un',0);

try 
    do_analysis = cell2mat(cellfun(@isempty,summary_table.Analyzed_ratio,'un',0));
catch
    do_analysis = isnan(summary_table.Analyzed_ratio);
end

keep_files = logical(do_analysis.*is_vClamp);
if any(~do_analysis)
    temp_disp = filenames(~do_analysis);
    disp(temp_disp(:));
    disp('found analysis files')
end


cellfun(@(x) int_pyr_ratio(x,false),filenames(keep_files),'un',0);

%%
clearvars -except PV_* SOM_*; clc;

% set plot params
area_names = {'LM','AL','PM','AM'};
setMin = true;
earlyLate_bound = 150;
cell_name = 'Pyr';
if strcmp(cell_name,'PV')
    int_color = [0.3 0.5 0.8];
elseif strcmp(cell_name,'SOM')
    int_color = [0.3 .8 0.4];
else
    int_color = [0.7 0.7 0.7];
end
cmap = brewermap(4,'Set2');
cmap_cell = mat2cell(cmap,[1 1 1 1],3);

% load data
excel_path = 'grouped_HVA_inhibition_fileList.xlsx';
load_path = 'Z:\home\jen\Analysis\HVA Recording\Int_pyr_ratio\';
summary_table = readtable(excel_path);
type_idx = find(sum(cell2mat(cellfun(@(x) contains(x,cell_name),[summary_table.HS1 summary_table.HS2],'un',0)),2));

% group data across expts 
EPSC_grouped = cell(4,1);
param_grouped =  cell(4,1);
collapse_minCurrent = [];
for area_i = 1:4
    idx = find(cell2mat(cellfun(@(x) contains(x,area_names{area_i}),summary_table.Area,'un',0)));
    idx = intersect(idx,type_idx);
       
    for file_i = 1:size(idx,1)
        load([load_path,summary_table.Date{idx(file_i)},summary_table.HFS{idx(file_i)},'.mat']);
        EPSC_grouped{area_i} = [EPSC_grouped{area_i} EPSC];
        param_grouped{area_i} = [param_grouped{area_i} param];
    end
    time_toPeak{area_i} = vertcat(EPSC_grouped{area_i}.timeToPeak)/20;
    int_pyr_ratio_grouped{area_i} = [EPSC_grouped{area_i}.Int_Pyr_ratio];
    filenames_grouped{area_i} = {param_grouped{area_i}.filename};
    
    laser_grouped{area_i} = [param_grouped{area_i}.laser_V];
    laser_grouped{area_i}(laser_grouped{area_i}>5) = 5;
    nCellsByArea{area_i} = numel(filenames_grouped{area_i});
    
    grouped_trace{area_i} = cat(3,EPSC_grouped{area_i}.meanTrace);
    denominator = min(grouped_trace{area_i}(:,2,:));
    temptrace = grouped_trace{area_i}./abs(denominator);
    norm_traces{area_i} = squeeze(temptrace);
    
    collapse_minCurrent = [collapse_minCurrent; vertcat(EPSC_grouped{area_i}.minCurrent)];
    
    temp_cat = vertcat(EPSC_grouped{area_i}.meanTraceRise);
    temp_cat(temp_cat>1) = temp_cat(temp_cat>1)/20000; % some expt have different timing type for some reason
    pyr_rise_time{area_i} = temp_cat(:,2)*1000;
    IN_rise_time{area_i} = temp_cat(:,1)*1000;
    
    temp_cat = vertcat(EPSC_grouped{area_i}.meanTraceLatency);
    pyr_latency{area_i} = temp_cat(:,2)./20;
    IN_latency{area_i} = temp_cat(:,1)./20;
    
    temp_cat = vertcat(EPSC_grouped{area_i}.jitter);
try
    if strcmp(cell_name,'SOM')
        IN_jitter{area_i} = temp_cat(:,1)*1000*.9;
    else
        IN_jitter{area_i} = temp_cat(:,1)*1000;
    end
end
    pyr_jitter{area_i} = temp_cat(:,2)*1000;
end
minCurrent = {collapse_minCurrent(:,1) collapse_minCurrent(:,2)};
timeVector = make_time(grouped_trace{area_i},20000,1);
temp = cellfun(@(x,y) repmat(x,y,1),area_names,nCellsByArea,'un',0);

traces_lat = cat(3,grouped_trace{1},grouped_trace{2});
traces_med = cat(3,grouped_trace{3},grouped_trace{4});

traces_pos = cat(3,grouped_trace{1},grouped_trace{3});
traces_ant = cat(3,grouped_trace{2},grouped_trace{4});

% calculate early and late EPSCs
for area_i = 1:4
    for file_i = 1:numel(EPSC_grouped{area_i})
        [temp_denom,temp_idx] = min(EPSC_grouped{area_i}(file_i).meanTrace(30:250,:,:));
        temp_idx = temp_idx + 30;
        if setMin
            temp_idx = EPSC_grouped{area_i}(file_i).earlyLatePoint;
        end
            temp_lat = EPSC_grouped{area_i}(file_i).meanTraceLatency;
            temp_trace_norm{area_i}(file_i,:,1) = rescale(-cumsum(EPSC_grouped{area_i}(file_i).meanTrace(temp_lat(1):temp_lat(1)+200,1)),0,1);
            temp_trace_norm{area_i}(file_i,:,2) = rescale(-cumsum(EPSC_grouped{area_i}(file_i).meanTrace(temp_lat(2):temp_lat(2)+200,2)),0,1);

            temp_trace{area_i}(file_i,:,:) = EPSC_grouped{area_i}(file_i).meanTrace;
        
        IN_early_EPSC{area_i}(file_i) = sum(temp_trace{area_i}(file_i,1:temp_idx(1),1));
        IN_late_EPSC{area_i}(file_i) = sum(temp_trace{area_i}(file_i,temp_idx(1):end,1));
        
        pyr_early_EPSC{area_i}(file_i) = sum(temp_trace{area_i}(file_i,1:temp_idx(2),2));
        pyr_late_EPSC{area_i}(file_i) = sum(temp_trace{area_i}(file_i,temp_idx(2):end,2));
        
        IN_early_EPSC_norm{area_i}(file_i) = (temp_trace_norm{area_i}(file_i,temp_idx(1),1));
        IN_late_EPSC_norm{area_i}(file_i) = 1-(temp_trace_norm{area_i}(file_i,temp_idx(1),1));
        
        pyr_early_EPSC_norm{area_i}(file_i) = (temp_trace_norm{area_i}(file_i,temp_idx(2),2));
        pyr_late_EPSC_norm{area_i}(file_i) = 1-(temp_trace_norm{area_i}(file_i,temp_idx(2),2));
        
        half_point_IN{area_i}(file_i) = find(temp_trace_norm{area_i}(file_i,:,1)>0.5,1,'first')./20;
        half_point_pyr{area_i}(file_i) = find(temp_trace_norm{area_i}(file_i,:,2)>0.5,1,'first')./20;
    end
end


nCellsByArea = cellfun(@numel,IN_early_EPSC,'un',0);

temp = cellfun(@(x,y) repmat(x,y,1),area_names,nCellsByArea,'un',0);
area_ids = vertcat(temp{:});

%% Figure 5C-D by cell type

fs = 10;
temp = cellfun(@(x,y) repmat(x,y,1),area_names,nCellsByArea,'un',0);
area_ids = vertcat(temp{:});

figure;

% onset latency
subplot(1,2,1);
if strcmp(cell_name,'Pyr')
    fast_errbar(1:4,pyr_latency,1,'continuous',false);
else
    fast_errbar(1:4,IN_latency,1,'continuous',false);
end
fix_axes(gcf,fs,'Area','Latency (ms)');
xticks(1:4);xticklabels(area_names);xlim([0 5]);
axis square; ylim([0 8]);

% 20-80 rise time
subplot(1,2,2);
if strcmp(cell_name,'Pyr')
    fast_errbar(1:4,pyr_rise_time,1,'continuous',false);
else
    fast_errbar(1:4,IN_rise_time,1,'continuous',false);
end
fix_axes(gcf,fs,'Area','20-80 rise time (ms)');
xticks(1:4);xticklabels(area_names);xlim([0 5]);
ylim([0 2]);
axis square;

%% Figure 5E-F

figure; hold on;
subplot(1,4,1);
plotSpread(int_pyr_ratio_grouped,'binWidth',1,'distributionMarkers','o','distributionColors',int_color);
fast_errbar(1:4,int_pyr_ratio_grouped,2,'continuous',false);
axis square;
xticks(1:4);
xticklabels(area_names');
xlim([0 5]);
fix_axes(gcf,fs,'Area','IN:Pyr');
title(cell_name);

%% Figure 6C by cell type
figure;
timeVector = make_time(temp_trace_norm{1},20,2);
fs = 10;
subplot(1,4,1); hold on;
for area_i = 1:4
if strcmp(cell_name,'Pyr')    
    plot(timeVector,mean(temp_trace_norm{area_i}(:,:,2)),'Color',cmap(area_i,:));
else
    plot(timeVector,mean(temp_trace_norm{area_i}(:,:,1)),'Color',cmap(area_i,:));
end
end
title([cell_name,' mean trace - norm.']);
axis tight;fix_axes(gcf,fs,'time (ms)','cumulative charge');
axis square

%% Figure 6D by cell type

all_hvas = {'LM','AL','PM','AM'};
fs = 10;
figure;
subplot(1,4,1)
if strcmp(cell_name,'Pyr')
    plotSpread(half_point_pyr,'distributionMarkers','o','binWidth',1,'distributionColors',int_color);
    fast_errbar(1:4,half_point_pyr,2,'continuous',false);
else
    plotSpread(half_point_IN,'distributionMarkers','o','binWidth',1,'distributionColors',int_color);
    fast_errbar(1:4,half_point_IN,2,'continuous',false);
end
axis square; fix_axes(gcf,fs,'Area','half point (s)');
title([cell_name,' half point']);
xticks(1:4);xticklabels(all_hvas); 
ylim([0 10]);

% stats 
% [p_half_IN,~,stats_result_half_IN] = anova1(cell2mat(half_point_IN)',area_ids,'off');
% multcompare(stats_result_half_IN)
% 
% [p_half_pyr,~,stats_result_half_pyr] = anova1(cell2mat(half_point_pyr)',area_ids,'off');
% multcompare(stats_result_half_pyr)

%% Figure 6E by cell type 

doNorm= true;
normSelf = true;
medLat = true;
if strcmp(cell_name,'Pyr')
    cell_type = 2;
else
    cell_type = 1;
end
    
if medLat
    traces_1 = traces_lat;
    traces_2 = traces_med;
    titles = [{'lateral'},{'medial'}];
else
    traces_1 = traces_ant;
    traces_2 = traces_pos;
    titles = [{'anterior'},{'posterior'}];
end
timeVector = make_time(traces_lat,20000,1);

figure;hold on;
if doNorm
    if normSelf
        fast_errbar(timeVector,squeeze(traces_1(:,cell_type,:)./abs(min(traces_1(30:120,cell_type,:)))),2,'shaded',true,'color',(int_color-0.3));
        fast_errbar(timeVector,squeeze(traces_2(:,cell_type,:)./abs(min(traces_2(30:120,cell_type,:)))),2,'shaded',true,'color',int_color);
    else
        fast_errbar(timeVector,squeeze(traces_1(:,cell_type,:)./abs(min(traces_1(30:120,2,:)))),2,'shaded',true,'color',int_color-0.3)
        fast_errbar(timeVector,squeeze(traces_2(:,cell_type,:)./abs(min(traces_2(30:120,2,:)))),2,'shaded',true,'color',(int_color))
    end
else
    fast_errbar(timeVector,squeeze(traces_1(:,cell_type,:)),2,'shaded',true,'color',(int_color-0.3))
    fast_errbar(timeVector,squeeze(traces_2(:,cell_type,:)),2,'shaded',true,'color',(int_color))
end
legend([titles(2);{''};titles(1);{''}]);
axis tight; axis square;
fix_axes(gcf,10,'Time (s)','amplitude');
title(cell_name);

%% Figure 6 - supplemental figure 1

for area_i = 1:4
    
    tempplot = vertcat(EPSC_grouped{area_i}.minCurrent);
    pyr_min{area_i} = tempplot(:,2);
end

figure;
% plotSpread(pyr_min,'distributionMarkers','o','binWidth',1,'distributionColors',cmap);
fast_errbar(1:4,pyr_min,1,'cells_as_x',true,'marker_size',15);
xticks(1:4);xticklabels(area_names); ylim([-1000 0]);
xlim([0 5]);
fix_axes(gcf,40,'Area','pA');
