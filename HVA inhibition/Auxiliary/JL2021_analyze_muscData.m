
function JL2021_analyze_muscData(input_datename,sheet_i,doPlot,your_save_path)
%% experiment metadata
if sheet_i == 1
    int_color = [0.3 .8 0.4];
else
    int_color = [0.3 0.5 0.8];
end

excel_filePath = 'muscimol_flowIn_PV_SOM.xlsx';
summary_table = readtable(excel_filePath,'Sheet',sheet_i);

file_idx = find(cell2mat(cellfun(@(x) strcmp(input_datename,x),summary_table.Date,'un',0)));
temp_date = strsplit(summary_table.Date{file_idx},'_');
filenames = strsplit(summary_table.Files{file_idx},',');
nFiles = numel(filenames);

if iscell(summary_table.LastSweep(file_idx))
    if isempty(summary_table.LastSweep{file_idx})
        lastSweeps = [];
    else
        lastSweeps = str2num(summary_table.LastSweep{file_idx});
    end
else
    lastSweeps = [];
end

if nFiles == 1 % all 3 conditions in one file
    filename = [temp_date{1},summary_table.Files{file_idx}];
    temp_dat = abfload([filename,'.abf']);
    muscIn_sweep = 10;
    muscOut_sweep = 30;
    baseline_sweeps = 1:10;
    drug_sweeps = 20:30;
elseif nFiles == 2 % + control and muscimol only
    [temp_input,temp_dat] = concat_sessions(summary_table.Date{file_idx},summary_table.Files{file_idx},lastSweeps);
    muscIn_sweep = temp_input.nTrials(1);
    muscOut_sweep = sum(temp_input.nTrials);
    baseline_sweeps = muscIn_sweep-10:muscIn_sweep;
    drug_sweeps = muscOut_sweep-15:muscOut_sweep;
elseif nFiles == 3 % control, flowin, +muscimol
    [temp_input,temp_dat] = concat_sessions(summary_table.Date{file_idx},summary_table.Files{file_idx},lastSweeps);
    muscIn_sweep = temp_input.nTrials(1);
    muscOut_sweep = sum(temp_input.nTrials);
    baseline_sweeps = muscIn_sweep-10:muscIn_sweep;
    drug_sweeps = muscOut_sweep-15:muscOut_sweep;
end

if contains(input_datename,'21121_3') % this day is weird timing
    baseline_sweeps = muscIn_sweep-5:muscIn_sweep;
    drug_sweeps = muscIn_sweep:muscIn_sweep+5;
end
    
nMins = 1:4:size(temp_dat,3);

if strcmp(summary_table.HS1{file_idx},'Pyr')
    IN_dim = 3;
    pyr_dim = 1;
else
    IN_dim = 1;
    pyr_dim = 3;
end


%% get laser on 
threshold = 0.3;
baseline = 0.05;
signal = 0.02;
sampleFreq = 20000;

norm_laser = squeeze(temp_dat(:,5,:))/5;
laser_on = norm_laser>threshold;
stimOn_idx = sum(cumprod(~laser_on));

% clean data
dat.IN = cleanData(temp_dat,IN_dim,stimOn_idx,'afterStimS',signal);
dat.Pyr = cleanData(temp_dat,pyr_dim,stimOn_idx,'afterStimS',signal);

%% show membrane r by time
threshold = 0.3;
baseline = 0.02;
signal = 0.03;
sampleFreq = 20000;

norm_laser = rescale(squeeze(abs(temp_dat(:,2,:))),0,1);
laser_on = norm_laser>threshold;
stimOn_idx = sum(cumprod(~laser_on));

% clean data
temp_IN_Rm = cleanData(temp_dat,IN_dim,stimOn_idx,'baselineS',baseline,'afterStimS',signal);
temp_Pyr_Rm = cleanData(temp_dat,pyr_dim,stimOn_idx,'baselineS',baseline,'afterStimS',signal);

dat.Pyr.Rm = (-.005./(mean(temp_Pyr_Rm.pts(end-100:end,:))*1E-12))*1E-6; % V/I = R
dat.IN.Rm = -.005./(mean(temp_IN_Rm.pts(end-100:end,:))*1E-12)*1E-6; % V/I = R

dat.Pyr.Rs = (-.005./(min(temp_Pyr_Rm.pts(1:10,:))*1E-12))*1E-6; % V/I = R
dat.IN.Rs = -.005./(min(temp_IN_Rm.pts(1:10,:))*1E-12)*1E-6; % V/I = R

% avg traces

dat.Pyr.meanTrace_baseline = mean(dat.Pyr.pts(:,baseline_sweeps),2);
[~,pyr_dI_min_idx] = min(dat.Pyr.meanTrace_baseline(1:100));
dat.Pyr.minCurrent_baseline = mean(dat.Pyr.meanTrace_baseline(pyr_dI_min_idx-5:pyr_dI_min_idx+5));
dat.Pyr.meanTraceLatency_baseline = (30+find(smooth(dat.Pyr.meanTrace_baseline(30:end))<dat.Pyr.minCurrent_baseline*0.2,1,'first'))/20000;

dat.IN.meanTrace_baseline = mean(dat.IN.pts(:,baseline_sweeps),2);
[~,IN_dI_min_idx] = min(dat.IN.meanTrace_baseline(1:100));
dat.IN.minCurrent_baseline = mean(dat.IN.meanTrace_baseline(IN_dI_min_idx-5:IN_dI_min_idx+5));
dat.IN.meanTraceLatency_baseline = (30+find(smooth(dat.IN.meanTrace_baseline(30:end))<dat.IN.minCurrent_baseline*0.2,1,'first'))/20000;

dat.Pyr.meanTrace_muscimol = mean(dat.Pyr.pts(:,drug_sweeps),2);
dat.Pyr.minCurrent_muscimol = mean(dat.Pyr.meanTrace_muscimol(pyr_dI_min_idx-5:pyr_dI_min_idx+5));
dat.Pyr.meanTraceLatency_muscimol = (30+find(smooth(dat.Pyr.meanTrace_muscimol(30:end))<dat.Pyr.minCurrent_muscimol*0.2,1,'first'))/20000;

dat.IN.meanTrace_muscimol = mean(dat.IN.pts(:,drug_sweeps),2);
dat.IN.minCurrent_muscimol = mean(dat.IN.meanTrace_muscimol(IN_dI_min_idx-5:IN_dI_min_idx+5));
dat.IN.meanTraceLatency_muscimol = (30+find(smooth(dat.IN.meanTrace_muscimol(30:end))<dat.IN.minCurrent_muscimol*0.2,1,'first'))/20000;

dat.Int_Pyr_ratio_baseline = dat.IN.minCurrent_baseline/dat.Pyr.minCurrent_baseline;
dat.Int_Pyr_ratio_muscimol = dat.IN.minCurrent_muscimol/dat.Pyr.minCurrent_muscimol;

% trial by trial
crop_dat_pyr = dat.Pyr.pts(pyr_dI_min_idx-5:pyr_dI_min_idx+5,:);
crop_dat_IN = dat.IN.pts(IN_dI_min_idx-5:IN_dI_min_idx+5,:);

for sweep_i = 1:size(temp_dat,3)
    dat.Pyr.tMinCurrent(sweep_i) = mean(crop_dat_pyr(:,sweep_i));
    dat.IN.tMinCurrent(sweep_i) = mean(crop_dat_IN(:,sweep_i));
end

%% plots
if doPlot
figure;
subplot(2,1,1); hold on;
plot(dat.Pyr.Rm,'k'); 
plot(dat.IN.Rm,'Color',int_color);
xticks(nMins);xticklabels(1:numel(nMins));
fix_axes(gcf,20,'time (min)','Rm'); axis tight
vline(muscIn_sweep,'k--');vline(muscOut_sweep,'k--');

subplot(2,1,2); hold on;
plot(dat.Pyr.Rs,'k'); 
plot(dat.IN.Rs,'Color',int_color);
xticks(nMins);xticklabels(1:numel(nMins));
fix_axes(gcf,20,'time (min)','Rs'); axis tight; ylim([0 50]);
vline(muscIn_sweep,'k--');vline(muscOut_sweep,'k--');
hline(30,'k--');

figure;
subplot(2,1,1);hold on;
plot(dat.Pyr.tMinCurrent,'ko');
plot(dat.IN.tMinCurrent,'o','Color',int_color);
xticks(nMins);xticklabels(1:numel(nMins));
vline(muscIn_sweep,'k--'); vline(muscOut_sweep,'k--');
fix_axes(gcf,10,'time (min)','pA');

subplot(2,1,2);hold on;
plot(dat.Pyr.tMinCurrent./dat.Pyr.tMinCurrent(1),'ko');
plot(dat.IN.tMinCurrent./dat.IN.tMinCurrent(1),'o','Color',int_color);
nMins = 1:4:size(temp_dat,3);
xticks(nMins);xticklabels(1:numel(nMins));
ylim([0 4]);
vline(muscIn_sweep,'k--'); vline(muscOut_sweep,'k--');

fix_axes(gcf,10,'time (min)','trial pA/first trial pA');

timeVector = make_time(dat.IN.pts,sampleFreq/1000,1);

figure;
subplot(1,2,1);hold on;
plot(timeVector,dat.Pyr.meanTrace_baseline,'k');
plot(timeVector,dat.Pyr.meanTrace_muscimol,'Color',[0.7 0.7 0.7]);
fix_axes(gcf,15,'time (ms)','pA'); axis square;
title('Pyr baseline vs muscimol')

subplot(1,2,2); hold on;
plot(timeVector,dat.IN.meanTrace_baseline,'Color',int_color);
plot(timeVector,dat.IN.meanTrace_muscimol,'Color',int_color+0.1);
fix_axes(gcf,15,'time (ms)','pA'); axis square;
title('IN baseline vs muscimol')


figure;
subplot(1,2,1);hold on;
plot([1 2],1000*[dat.Pyr.meanTraceLatency_baseline dat.Pyr.meanTraceLatency_muscimol],'k');
plot([1 2],1000*[dat.IN.meanTraceLatency_baseline dat.IN.meanTraceLatency_muscimol],'Color',int_color);
ylim([0 5]);xticks(1:2);xticklabels({'baseline';'drug'});
xlim([0.5 2.5]);
fix_axes(gcf,15,'Condition','Latency diff'); axis square;

subplot(1,2,2); hold on;
plot([1 2],[dat.Int_Pyr_ratio_baseline dat.Int_Pyr_ratio_muscimol],'k');
% ylim([-0.2 0.2]);
xticks(1:2);xticklabels({'baseline';'drug'});
xlim([0.5 2.5]);ylim([0 5])
fix_axes(gcf,15,'Condition','IN:Pyr ratio');axis square;
end
%%

save_path = [your_save_path,summary_table.Date{file_idx}];
save(save_path,'dat');
notebook_path = 'muscimol_flowIn_PV_SOM.xlsx';
try
    xlswrite(notebook_path,...
    {save_path},sheet_i,...
    ['K' num2str(file_idx + 1)]);
    xlswrite(notebook_path,...
    {'X'},sheet_i,['J' num2str(file_idx + 1)]);
catch
    input('close the excel file. press enter when done.');
    xlswrite(notebook_path,...
    {save_path},sheet_i,...
    ['K' num2str(file_idx + 1)]);
    xlswrite(notebook_path,...
    {'X'},sheet_i,['J' num2str(file_idx + 1)]);
end
end