%% code for recreating figures 6 and 7 of Li2021
your_read_path = ['Z:\home\jen\Notebook files\HVA recording\']; % change for your path
your_save_path = ['Z:\home\jen\Analysis\HVA Recording\Int_pyr_ratio\']; % change for your path
%%
JL2021_load_intPyr_ratio(your_read_path,your_save_path);
PV = load([your_save_path,'PV_grouped.mat']);
SOM = load([your_save_path,'SOM_grouped.mat']);
Pyr = load([your_save_path,'Pyr_grouped.mat']);


cell_names = {'PV','SOM','Pyr'};
area_names = {'LM','AL','PM','AM'};
int_color = [0.3 0.5 0.8;
    0.3 .8 0.4;
    0.7 0.7 0.7];
cmap = brewermap(4,'Set2');
cmap_cell = mat2cell(cmap,[1 1 1 1],3);
disp('analyzed and loaded');

%% Figure 5H
doStats = true;
fs = 10;
figure; hold on;
subplot(1,4,1);
plotSpread(cellfun(@(x) log10(x),PV.int_pyr_ratio_grouped,'un',0),'binWidth',1,'distributionMarkers','o','distributionColors',int_color(1,:));
fast_errbar(1:4,cellfun(@(x) log10(x),PV.int_pyr_ratio_grouped,'un',0),2,'continuous',false,'stats',doStats);
axis square;
xticks(1:4);
xticklabels(area_names');
xlim([0 5]);
fix_axes(gcf,fs,'Area','PV:Pyr');
ylim([-1 2]);yticks(-1:0.5:2)
subplot(1,4,2);
plotSpread(cellfun(@(x) log10(x),SOM.int_pyr_ratio_grouped,'un',0),'binWidth',1,'distributionMarkers','o','distributionColors',int_color(2,:));
fast_errbar(1:4,cellfun(@(x) log10(x),SOM.int_pyr_ratio_grouped,'un',0),2,'continuous',false,'stats',doStats);
axis square;
xticks(1:4);
xticklabels(area_names');
xlim([0 5]);
ylim([-2 1]);yticks(-2:0.5:1);

fix_axes(gcf,fs,'Area','SOM:Pyr');

%% Figure 6C by cell type
figure;
timeVector = make_time(PV.cumul_trace_norm{1},20,2);
fs = 10;
for area_i = 1:4
    subplot(1,4,1); hold on;
    plot(timeVector,mean(squeeze(Pyr.cumul_trace_norm{area_i}(:,:,2))),'Color',cmap(area_i,:));
    title('Pyr mean trace - norm.');
    axis tight;fix_axes(gcf,fs,'time (ms)','cumulative charge');
    axis square
    
    subplot(1,4,2); hold on;
    plot(timeVector,mean(squeeze(PV.cumul_trace_norm{area_i}(:,:,1))),'Color',cmap(area_i,:));
    title('PV mean trace - norm.');
    axis tight;fix_axes(gcf,fs,'time (ms)','cumulative charge');
    axis square
    
    subplot(1,4,3); hold on;
    plot(timeVector,mean(squeeze(SOM.cumul_trace_norm{area_i}(:,:,1))),'Color',cmap(area_i,:));
    title('SOM mean trace - norm.');
    axis tight;fix_axes(gcf,fs,'time (ms)','cumulative charge');
    axis square
end

%% Figure 6D by cell type

all_hvas = {'LM','AL','PM','AM'};
fs = 10;
figure;
subplot(1,4,1)
plotSpread(Pyr.half_point_pyr,'distributionMarkers','o','binWidth',1,'distributionColors',int_color(3,:));
fast_errbar(1:4,Pyr.half_point_pyr,2,'continuous',false);
axis square; fix_axes(gcf,fs,'Area','half point (ms)');
title(['Pyr half point']);
xticks(1:4);xticklabels(all_hvas); 
ylim([0 10]);

subplot(1,4,2);
plotSpread(PV.half_point_IN,'distributionMarkers','o','binWidth',1,'distributionColors',int_color(1,:));
fast_errbar(1:4,PV.half_point_IN,2,'continuous',false);
axis square; fix_axes(gcf,fs,'Area','half point (ms)');
title(['PV half point']);
xticks(1:4);xticklabels(all_hvas); 
ylim([0 10]);

subplot(1,4,3);
plotSpread(SOM.half_point_IN,'distributionMarkers','o','binWidth',1,'distributionColors',int_color(2,:));
fast_errbar(1:4,SOM.half_point_IN,2,'continuous',false);
axis square; fix_axes(gcf,fs,'Area','half point (ms)');
title(['SOM half point']);
xticks(1:4);xticklabels(all_hvas); 
ylim([0 10]);

%% Figure 6E by cell type 
doNorm= true;
normSelf = true;
medLat = true;
timeVector = make_time(Pyr.traces_lat,20000,1);


plotStruct = [PV SOM Pyr];
cell_types = [1 1 2];

figure;
for i_type = 1:3
    cell_type = cell_types(i_type);
    subplot(1,3,i_type); hold on;
if medLat
    traces_1 = plotStruct(i_type).traces_lat;
    traces_2 = plotStruct(i_type).traces_med;
    titles = [{'lateral'},{'medial'}];
else
    traces_1 = plotStruct(i_type).traces_ant;
    traces_2 = plotStruct(i_type).traces_pos;
    titles = [{'anterior'},{'posterior'}];
end

if doNorm
    if normSelf
        fast_errbar(timeVector,squeeze(traces_1(:,cell_type,:)./abs(min(traces_1(30:120,cell_type,:)))),2,'shaded',true,'color',(int_color(i_type,:)-0.3));
        fast_errbar(timeVector,squeeze(traces_2(:,cell_type,:)./abs(min(traces_2(30:120,cell_type,:)))),2,'shaded',true,'color',int_color(i_type,:));
    else
        fast_errbar(timeVector,squeeze(traces_1(:,cell_type,:)./abs(min(traces_1(30:120,2,:)))),2,'shaded',true,'color',int_color(i_type,:)-0.3)
        fast_errbar(timeVector,squeeze(traces_2(:,cell_type,:)./abs(min(traces_2(30:120,2,:)))),2,'shaded',true,'color',(int_color(i_type,:)))
    end
else
    fast_errbar(timeVector,squeeze(traces_1(:,cell_type,:)),2,'shaded',true,'color',(int_color(i_type,:)-0.3))
    fast_errbar(timeVector,squeeze(traces_2(:,cell_type,:)),2,'shaded',true,'color',(int_color(i_type,:)))
end
legend([titles(2);{''};titles(1);{''}]);
axis tight; axis square;
fix_axes(gcf,10,'Time (s)','amplitude');
title(cell_names{i_type});
end
%% Supplementary figure 4C
load([your_save_path,'interneurons_multipower.mat'])

inpow_pop = [];
cell_type = {'all_pv','all_som','py_l23'};
for i_ex = 1:numel(dat)
    
    
    % both recording channels should have valid Vclamp data. Check this,
    % and then (optionally) exclude files with very small peak epscs.
    assert(all(dat{i_ex}.info.HS_is_valid_Vclamp), 'Error: vclamp not defined for at least one channel')

    %
    % analyze the first pulse amplitude across powers
    %
    %%%%%%%
    inpow_pop.dat{i_ex}.p1_amp = {};
    inpow_pop.dat{i_ex}.laser_V = [];
    condnames = fieldnames(dat{i_ex}.expt);
    ex_laser_powers = cellfun(@(x) dat{i_ex}.expt.(x).tdict(1), condnames);
    l_ex_good_powers = true(size(ex_laser_powers));
    ch_amps = {};
    for i_ch = 1:2
        amps = cellfun(@(x) dat{i_ex}.expt.(x).stats.EPSCamp{i_ch}(1,:,:), condnames, 'uniformoutput', false);
        % hack to remove some sweeps that were deleted in the .abf file
        % post-hoc. Greedy update of exp-wide list
        l_ch_good_powers = cellfun(@(x) ~isempty(x), amps);
        l_ex_good_powers = l_ex_good_powers & l_ch_good_powers;
        ch_amps{i_ch} = amps;
    end
    for i_ch = 1:2
        inpow_pop.dat{i_ex}.p1_amp{i_ch} = cat(1, ch_amps{i_ch}{l_ex_good_powers});
    end
    inpow_pop.dat{i_ex}.laser_V = ex_laser_powers(l_ex_good_powers);
    
    above_10pa = cellfun(@(x) x>10, inpow_pop.dat{i_ex}.p1_amp, 'uniformoutput', false);
    l_gt_10pA = above_10pa{1} & above_10pa{2};
    inpow_pop.dat{i_ex}.p1_amp = cellfun(@(x) x(l_gt_10pA), inpow_pop.dat{i_ex}.p1_amp, 'uniformoutput', false);
    inpow_pop.dat{i_ex}.laser_V = inpow_pop.dat{i_ex}.laser_V(l_gt_10pA);
    
    
    % grab some meta data
    inpow_pop.dat{i_ex}.hva = dat{i_ex}.info.brainArea;
    inpow_pop.dat{i_ex}.cell_type = dat{i_ex}.info.cellType;
end

PLOT_TYPE = 'ratio'; % could be 'ratio', or 'scatter'
COLOR_BY = 'in';  % could be 'hva', or 'in'

NORM_TO_PY_MAX = false;
COMBINE_HVAS = true;
COMBINE_INS = true;

assert(strcmp(EXPTTYPE, 'IN_powers'), 'ERROR: not the correct type of experiment for multipower analysis')
assert(COMBINE_INS, 'error: combine interneurons needs to be true, otherwise things may break')

figure, hold on,

for i_ex = 1:numel(dat)
    
    tmp_dat = {[], []};
    for i_ch = 1:2
        
        ex_brain_area = inpow_pop.dat{i_ex}.hva;
        if COMBINE_HVAS
            if regexpi(ex_brain_area, 'pm|am')
                ex_brain_area = 'med';
                hva_clr = 'b';
            elseif regexpi(ex_brain_area, 'lm|al')
                ex_brain_area = 'lat';
                hva_clr = 'r';
            else
                error('did not identify area')
            end
        end
        
        temp_cell_type = inpow_pop.dat{i_ex}.cell_type{i_ch};
        if ~strcmpi(temp_cell_type, 'py_l23')
            if COMBINE_INS
                if regexpi(temp_cell_type, 'LTSIN|SOMCRE')
                    cell_type{i_ex} = 'all_som';
                    cell_clr = int_color(2,:);
                elseif regexpi(temp_cell_type, 'FS|PVCRE')
                    cell_type{i_ex} = 'all_pv';
                    cell_clr = int_color(1,:);
                else
                    error('did not identify cell type')
                end
            end
            plt_idx = 2; % interneurons are on Y
        else
            plt_idx = 1;  % PY cell on X
        end
       
        tmp_dat{plt_idx} = inpow_pop.dat{i_ex}.p1_amp{i_ch};
    end
    
    if NORM_TO_PY_MAX
        maxval = tmp_dat{1}(end);
        tmp_dat = cellfun(@(x) x./maxval, tmp_dat, 'uniformoutput', false);
    end
    
    
    switch COLOR_BY
        case 'in'
            plt_clr = cell_clr;
        case 'hva'
            plt_clr = hva_clr;
    end
    
    % to plot the raw values
    if strcmpi(PLOT_TYPE, 'scatter')
        plot(tmp_dat{1}, tmp_dat{2}, '-o', 'color', plt_clr, 'linewidth', 2)
    elseif strcmpi(PLOT_TYPE, 'ratio')
        % to plot the ratio values
        subplot(1,2,1); hold on;
        ratio_vals{i_ex} = tmp_dat{2} ./ tmp_dat{1};
        plot(inpow_pop.dat{i_ex}.laser_V, ratio_vals{i_ex}, '-o', 'color', plt_clr, 'linewidth', 2)
        set(gca, 'yscale', 'log')
        ylabel('IN : PY ratio', 'fontsize', 14)
        xlabel('Laser Voltage', 'fontsize', 14)
        laser_powerV{i_ex} = inpow_pop.dat{i_ex}.laser_V;
        
    end
end
fix_axes(gcf,15,'','')
ylabel('IN : PY ratio', 'fontsize', 14)
xlabel('Laser Voltage', 'fontsize', 14)
axis square
%% stats for ratio vals
pv_idx = cell2mat(cellfun(@(x) strcmp(x,'all_pv'),cell_type,'un',0));
nPV = sum(pv_idx);
nSOM = sum(~pv_idx);
[pv_p,~,pv_stats] = anovan(cell2mat(ratio_vals'),cell2mat(laser_powerV'));


%% Supplementary figure 4D

doStats = true;
figure;
tempPlot = cellfun(@(x) x(:,2),Pyr.minCurrent_byArea,'un',0);
fast_errbar(1:4,tempPlot,1,'marker_size',15,'stats',doStats);
xticks(1:4);xticklabels(area_names); ylim([-1000 0]);
xlim([0 5]);
fix_axes(gcf,40,'Area','pA');

