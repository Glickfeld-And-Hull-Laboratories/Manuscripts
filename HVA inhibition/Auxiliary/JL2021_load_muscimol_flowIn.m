function JL2021_load_muscimol_flowIn(your_read_path,your_save_path);
%% analyze all data and save mat files

% read excel file and analyze files that don't have X in "Analyzed" column
excel_filePath = 'muscimol_flowIn_PV_SOM.xlsx';

for sheet_i = 1:2
summary_table = readtable(excel_filePath,'Sheet',sheet_i);

try 
    do_analysis = logical(cell2mat(cellfun(@isempty,summary_table.Analyzed,'un',0)));
catch
    do_analysis = logical(isnan(summary_table.Analyzed));
end

if any(~do_analysis)
    temp_disp = summary_table.Date(~do_analysis);
    disp([temp_disp(:)]);
    disp('found analysis files')
end

cellfun(@(x) analyze_muscData(x,sheet_i,false),summary_table.Date(do_analysis),'un',0);
end
%% load data  
cell_names = {'SOM';'PV';'Pyr'};
excel_filePath = 'Z:\home\jen\Notebook files\HVA Recording\muscimol_flowIn_PV_SOM.xlsx';
for i_type = 1:3
    summary_table = [];
    if i_type == 3
        for temp_sheet = 1:2
        temp_table = readtable(excel_filePath,'Sheet',temp_sheet);
        summary_table = [summary_table; temp_table];
        end
    else
        summary_table = readtable(excel_filePath,'Sheet',i_type);
    end

idx = find(cell2mat(cellfun(@(x) strcmp(x,'X'),summary_table.Keep,'un',0)));
data_grouped = [];
for file_i = 1:size(idx,1)
    load([summary_table.SavePath{idx(file_i)},'.mat']);
    data_grouped = [data_grouped dat];
    
    
    temp_lat = floor(dat.Pyr.meanTraceLatency_baseline*20000);
    pyr_cumul_baseline(file_i,:) = rescale(-cumsum(dat.Pyr.meanTrace_baseline(temp_lat:temp_lat+200,1)),0,1);
    pyr_cumul_musc(file_i,:) = rescale(-cumsum(dat.Pyr.meanTrace_muscimol(temp_lat:temp_lat+200,1)),0,1);
    
    temp_lat = floor(dat.IN.meanTraceLatency_baseline*20000);
    IN_cumul_baseline(file_i,:) = rescale(-cumsum(dat.IN.meanTrace_baseline(temp_lat:temp_lat+200,1)),0,1);
    IN_cumul_musc(file_i,:) = rescale(-cumsum(dat.IN.meanTrace_muscimol(temp_lat:temp_lat+200,1)),0,1);
end
temp_cat = vertcat(data_grouped.Pyr);
pyr_latency_baseline = vertcat(temp_cat.meanTraceLatency_baseline);
pyr_latency_muscimol = vertcat(temp_cat.meanTraceLatency_muscimol);

temp_cat = vertcat(data_grouped.IN);
IN_latency_baseline = vertcat(temp_cat.meanTraceLatency_baseline);
IN_latency_muscimol = vertcat(temp_cat.meanTraceLatency_muscimol);

Int_pyr_ratio_baseline = vertcat(data_grouped.Int_Pyr_ratio_baseline);
Int_pyr_ratio_muscimol = vertcat(data_grouped.Int_Pyr_ratio_muscimol);

temp_struct = [data_grouped.Pyr];
pyr_baselineTrace = [temp_struct.meanTrace_baseline];
pyr_muscimolTrace = [temp_struct.meanTrace_muscimol];

temp_struct = [data_grouped.IN];
IN_baselineTrace = [temp_struct.meanTrace_baseline];
IN_muscimolTrace = [temp_struct.meanTrace_muscimol];


clear temp* dat* summary_table
save([your_save_path,cell_names{i_type},'_grouped.mat']);
end

disp('loaded')
