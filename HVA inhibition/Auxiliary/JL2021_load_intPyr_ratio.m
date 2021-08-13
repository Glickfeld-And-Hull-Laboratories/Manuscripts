function JL2021_load_intPyr_ratio(your_read_path,your_save_path)

addpath(your_save_path);
summary_table = readtable([your_read_path,'Interneuron_properties_ratio.xlsx']);
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

cellfun(@(x) JL2021_int_pyr_ratio(x,false,your_save_path),filenames(keep_files),'un',0);

%% 
% set plot params
area_names = {'LM','AL','PM','AM'};
cell_names = {'PV','SOM','Pyr'};
for i_type = 1:3
    cell_name = cell_names{i_type};
    % load data
    summary_table = readtable([your_read_path,'Interneuron_properties_ratio.xlsx']);
    type_idx = find(sum(cell2mat(cellfun(@(x) contains(x,cell_name),[summary_table.HS1 summary_table.HS2],'un',0)),2));

    % group data across expts 
    EPSC_grouped = cell(4,1);
    param_grouped =  cell(4,1);
    temp_minCurrent = [];
    for area_i = 1:4
        temp_idx = find(cell2mat(cellfun(@(x) contains(x,area_names{area_i}),summary_table.Area,'un',0)));
        temp_idx = intersect(temp_idx,type_idx);
        keep_idx = find(cell2mat(cellfun(@(x) contains(x,'X'),summary_table.KeepRatio,'un',0)));
        temp_idx = intersect(temp_idx,keep_idx);

        for file_i = 1:size(temp_idx,1)
            load([your_save_path,summary_table.Date{temp_idx(file_i)},summary_table.HFS{temp_idx(file_i)},'.mat']);
            EPSC_grouped{area_i} = [EPSC_grouped{area_i} EPSC];
            param_grouped{area_i} = [param_grouped{area_i} param];
        end
        int_pyr_ratio_grouped{area_i} = [EPSC_grouped{area_i}.Int_Pyr_ratio];
        filenames_grouped{area_i} = {param_grouped{area_i}.filename};

        laser_grouped{area_i} = [param_grouped{area_i}.laser_V];
        laser_grouped{area_i}(laser_grouped{area_i}>5) = 5;
        nCellsByArea{area_i} = numel(filenames_grouped{area_i});

        grouped_trace{area_i} = cat(3,EPSC_grouped{area_i}.meanTrace);  
        temp_minCurrent = [temp_minCurrent; vertcat(EPSC_grouped{area_i}.minCurrent)];

        minCurrent_byArea{area_i} = vertcat(EPSC_grouped{area_i}.minCurrent);
        temp_cat = vertcat(EPSC_grouped{area_i}.meanTraceRise);
        temp_cat(temp_cat>1) = temp_cat(temp_cat>1)/20000; % some expt have different timing sample freq for some reason
        pyr_rise_time{area_i} = temp_cat(:,2)*1000;
        IN_rise_time{area_i} = temp_cat(:,1)*1000;

        temp_cat = vertcat(EPSC_grouped{area_i}.meanTraceLatency);
        pyr_latency{area_i} = temp_cat(:,2)./20;
        IN_latency{area_i} = temp_cat(:,1)./20;
    
    end
    minCurrent = {temp_minCurrent(:,1) temp_minCurrent(:,2)};
    traces_lat = cat(3,grouped_trace{1},grouped_trace{2});
    traces_med = cat(3,grouped_trace{3},grouped_trace{4});

    traces_pos = cat(3,grouped_trace{1},grouped_trace{3});
    traces_ant = cat(3,grouped_trace{2},grouped_trace{4});

    % calculate early and late EPSCs
    for area_i = 1:4
        for file_i = 1:numel(EPSC_grouped{area_i})
            [~,temp_idx] = min(EPSC_grouped{area_i}(file_i).meanTrace(30:250,:,:));
            temp_idx = temp_idx + 30;
            temp_lat = EPSC_grouped{area_i}(file_i).meanTraceLatency;
            cumul_trace_norm{area_i}(file_i,:,1) = rescale(-cumsum(EPSC_grouped{area_i}(file_i).meanTrace(temp_lat(1):temp_lat(1)+200,1)),0,1);
            cumul_trace_norm{area_i}(file_i,:,2) = rescale(-cumsum(EPSC_grouped{area_i}(file_i).meanTrace(temp_lat(2):temp_lat(2)+200,2)),0,1);

            half_point_IN{area_i}(file_i) = find(cumul_trace_norm{area_i}(file_i,:,1)>0.5,1,'first')./20;
            half_point_pyr{area_i}(file_i) = find(cumul_trace_norm{area_i}(file_i,:,2)>0.5,1,'first')./20;
        end
    end


    nCellsByArea = cellfun(@numel,half_point_IN,'un',0);

    temp = cellfun(@(x,y) repmat(x,y,1),area_names,nCellsByArea,'un',0);
    area_ids = vertcat(temp{:});

    clear temp* EPSC_grouped

    save([your_save_path,cell_name,'_grouped']);
end

