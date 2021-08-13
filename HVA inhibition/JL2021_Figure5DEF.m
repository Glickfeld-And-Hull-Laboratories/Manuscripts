ccc;
your_read_path = 'Z:\home\jen\Notebook files\HVA Recording\';
your_save_path = 'Z:\home\jen\Analysis\HVA Recording\Reviewer expt\';
addpath(genpath(your_read_path));
% JL2021_load_muscimol_flowIn(your_read_path,your_save_path);
%%
PV = load([your_save_path,'PV_grouped']);
SOM = load([your_save_path,'SOM_grouped']);
Pyr = load([your_save_path,'Pyr_grouped']);
cell_names = {'PV';'SOM';'Pyr'};

int_color = [0.3 0.5 0.8;
    0.3 .7 0.2;
    0.7 0.7 0.7];

%% Figure 5D
timeVector = make_time(Pyr.pyr_baselineTrace,20,1);
plotStruct = [PV; SOM; Pyr];
figure;
subplot(1,4,1)
i_type = 3;
denom = min([plotStruct(i_type).pyr_baselineTrace]);
fast_errbar(timeVector,-[plotStruct(i_type).pyr_baselineTrace]./denom,2,'shaded',true,'Color',int_color(i_type,:));
fast_errbar(timeVector,-[plotStruct(i_type).pyr_muscimolTrace]./denom,2,'shaded',true,'Color',int_color(i_type,:)+0.15);
fast_errbar(timeVector,-([plotStruct(i_type).pyr_baselineTrace]-[plotStruct(i_type).pyr_muscimolTrace])./denom,2,'shaded',true,'Color',int_color(i_type,:)*0.3);
fix_axes(gcf,10,'time (ms)','pA');vline(4);
ylim([-Inf 0.04]); axis square;

for i_type = 1:2
    subplot(1,4,i_type+1)
    denom = min([plotStruct(i_type).IN_baselineTrace]);
    fast_errbar(timeVector,-[plotStruct(i_type).IN_baselineTrace]./denom,2,'shaded',true,'Color',int_color(i_type,:));
    fast_errbar(timeVector,-[plotStruct(i_type).IN_muscimolTrace]./denom,2,'shaded',true,'Color',int_color(i_type,:)+0.15);
    fast_errbar(timeVector,-([plotStruct(i_type).IN_baselineTrace]-[plotStruct(i_type).IN_muscimolTrace])./denom,2,'shaded',true,'Color',int_color(i_type,:)*0.3);
    fix_axes(gcf,10,'time (ms)','pA');vline(4);
    ylim([-Inf 0.04]); axis square;
end

%% figure 5E-F
plotStruct = [PV; SOM; Pyr];
figure;
subplot(1,4,1); hold on;
for i_type = 1:2
    latency_diff{i_type} = 1000*([plotStruct(i_type).IN_latency_muscimol]-[plotStruct(i_type).IN_latency_baseline]);
end
latency_diff{3} = 1000*([plotStruct(i_type).pyr_latency_muscimol]-[plotStruct(i_type).pyr_latency_baseline]);
plotSpread(latency_diff,'distributionMarkers','o','distributionColors',int_color)
fast_errbar(4,cell2mat(latency_diff'),1);
ylim([-1 1]); xticks(1:4);xticklabels([cell_names;{'mean'}]);
hline(0,'k--');
fix_axes(gcf,10,'cell type','delta latency (ms)');
axis square;

ylabel_str = '% change IN:Pyr';
ylims = [-100 100];
subplot(1,4,2);
for i_type = 1:2
in_pyr_ratio_diff = 100*([plotStruct(i_type).Int_pyr_ratio_muscimol]-[plotStruct(i_type).Int_pyr_ratio_baseline])./[plotStruct(i_type).Int_pyr_ratio_baseline];
plot(i_type,in_pyr_ratio_diff,'o','Color',int_color(i_type,:)); hold on;
fast_errbar(i_type,in_pyr_ratio_diff,1);
end
ylim(ylims); xticks(1:2);xticklabels({'PV';'SOM'});xlim([0.5 2.5]);hline(0,'k--');
fix_axes(gcf,10,'',ylabel_str);
axis square;

clear plotStruct