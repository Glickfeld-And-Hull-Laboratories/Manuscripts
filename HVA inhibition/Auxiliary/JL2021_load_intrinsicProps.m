function JL2021_load_intrinsicProps(your_save_path)

%% analyze all data and save mat files

% read excel file and analyze files that don't have X in "Analyzed" column
for type_i = 1:3
summary_table = readtable('grouped_HVA_intrinsicProps_fileList.xlsx','Sheet',type_i);
filenames = cellfun(@(x,y) [x,y],summary_table.Date,summary_table.DCSteps,'un',0);

try 
    do_analysis = logical(cell2mat(cellfun(@isempty,summary_table.Analyzed,'un',0)));
catch
    do_analysis = logical(isnan(summary_table.Analyzed));
end

if any(~do_analysis)
        temp_disp = filenames(~do_analysis);
    disp([temp_disp(:)]);
    disp('found analysis files')
end

cellfun(@(x) JL2021_intrinsic_props(x,false),filenames(do_analysis),'un',0);
end
%% load data grouped by cell type and plot 
area_names = {'LM','AL','PM','AM'};
cell_names = {'PV','SOM','Pyr'};
level_idx = [2 2 1];
for i_type = 1:3
% load data
summary_table = readtable('grouped_HVA_intrinsicProps_fileList.xlsx','Sheet',i_type);

%%
% pre-allocate 
AP_adapt = cell(1,4);
half_width = cell(1,4);
spikes_grouped = cell(1,4);
param_grouped = cell(1,4);

% load into grouped struct and count cells for each expt type
for area_i = 1:4
    idx = find(cell2mat(cellfun(@(x) contains(x,area_names{area_i}),summary_table.Area,'un',0)));
    for file_i = 1:size(idx,1)
        load([your_save_path,summary_table.Date{idx(file_i)},summary_table.DCSteps{idx(file_i)},'.mat']);        
        HS_i = summary_table.HS_i(idx(file_i));
        
        if iscell(HS_i)
            HS_i = str2num(HS_i{1});
        end
        for cell_i = HS_i
            spikes_grouped{area_i} = [spikes_grouped{area_i} spikes(cell_i)];
            param_grouped{area_i} = [param_grouped{area_i} param];

            % check if there are spikes 
            temp_idx = level_idx(i_type)+find(~isnan(spikes(cell_i).AP_adapt),1,'first');
            if ~isempty(temp_idx) && temp_idx <= numel(param.I_vals)
                AP_adapt{area_i} = [AP_adapt{area_i} spikes(cell_i).AP_adapt(temp_idx)];
                half_width{area_i} = [half_width{area_i} spikes(cell_i).half_width(temp_idx)];
            else
                AP_adapt{area_i} = [AP_adapt{area_i} NaN];
                half_width{area_i} = [half_width{area_i} NaN];
            end
        end
    end

    FR_grouped{area_i} = [spikes_grouped{area_i}.FR];
    tau_grouped{area_i} = 1000*[spikes_grouped{area_i}.tau_63];
    Rin_grouped{area_i} = [spikes_grouped{area_i}.Rin];
    nCellsByArea{area_i} = sum(~isnan(Rin_grouped{area_i}));
    filenames{area_i} = {param_grouped{area_i}.filename};
end

% calculate FI and Ih for matched bins
posBins = [0,100,200,300,400,500];
if i_type == 3
    negBins = [-680,-450,-200,-100,0];
else
    negBins = [-400,-300,-200,-100,0];
end
for area_i = 1:4
    for cell_i = 1:numel(spikes_grouped{area_i})
        Ih_slope{area_i}(cell_i) = max(1000*diff(spikes_grouped{area_i}(cell_i).Ih_Vm(1:3))./diff(param_grouped{area_i}(cell_i).I_vals(1:3)));
        temp_i = discretize(posBins,param_grouped{area_i}(cell_i).I_vals,'IncludedEdge','left');
        for i_inj = 1:numel(posBins)
            if sum(~isnan(temp_i(i_inj)))>0
                FI{area_i}(cell_i,i_inj) = nanmean(spikes_grouped{area_i}(cell_i).FR(temp_i(i_inj)));
            else
                FI{area_i}(cell_i,i_inj) = NaN;
            end 
        end
        temp_i = discretize(negBins,param_grouped{area_i}(cell_i).I_vals);
        for i_inj = 1:numel(negBins)
            if sum(~isnan(temp_i(i_inj))>0)
                tau_temp(cell_i,i_inj) = nanmean(spikes_grouped{area_i}(cell_i).tau_63(temp_i(i_inj)));
                Ih{area_i}(cell_i,i_inj) = -nanmean(spikes_grouped{area_i}(cell_i).Ih_Vm(temp_i(i_inj)));
            else
                tau_temp(cell_i,i_inj) = NaN;
                Ih{area_i}(cell_i,i_inj) = NaN; 
            end
        end
        tau{area_i} = 1000*nanmean(tau_temp(:,[2,3]),2);
        
    end
    clear tau_temp    
    nCellsByArea_FI{area_i} = sum(sum(~isnan(FI{area_i}),2)>0);
end
%% run stats

temp = cellfun(@(x,y) repmat(x,y,1),area_names,nCellsByArea,'un',0);
area_ids = vertcat(temp{:});
[p_Rin,~,stats_result_Rin] = anova1(cell2mat(Rin_grouped)',area_ids,'off');
[p_HW,~,stats_result_HW] = anova1(cell2mat(half_width)',area_ids,'off');
[p_APadapt,~,stats_result_APadapt] = anova1(cell2mat(AP_adapt)',area_ids,'off');
[p_tau,~,stats_result_tau] = anova1(cell2mat(tau'),area_ids,'off');
[p_Ih,~,stats_result_Ih] = anova1(cell2mat(Ih_slope),area_ids,'off');

temp_mat = cell2mat(FI');
area_label = cellfun(@(x,y) repmat(x,y,size(temp_mat,2)),[{1},{2},{3},{4}],nCellsByArea_FI,'un',0);
area_label = vertcat(area_label{:});
Icmd = repmat(1:numel(posBins),size(temp_mat,1),1);
[p_FI,~,stats_result_FI] = anovan(temp_mat(:),[{area_label(:)},{Icmd(:)}],'varnames',{'area','Icmd'},'model','interaction','display','off');

%% collapsed cell characteristics for generating Table 1

all_Rin_mean = nanmean(cell2mat(Rin_grouped));
all_Rin_std = nanstd(cell2mat(Rin_grouped));

all_HW_mean = nanmean(cell2mat(half_width));
all_HW_std = nanstd(cell2mat(half_width));

all_adapt_mean = nanmean(cell2mat(AP_adapt));
all_adapt_std = nanstd(cell2mat(AP_adapt));

all_tau_mean = nanmean(cell2mat(tau'));
all_tau_std = nanstd(cell2mat(tau'));

all_Ihslope_mean = nanmean(cell2mat(Ih_slope));
all_Ihslope_std = nanstd(cell2mat(Ih_slope));

clear temp* spikes* param *_i do_analysis
save([your_save_path,cell_names{i_type},'_grouped.mat']);
clearvars -except i_type your_save_path cell_names level_idx area_names
end