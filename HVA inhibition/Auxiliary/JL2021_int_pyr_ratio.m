function [EPSC,param] = JL2021_int_pyr_ratio(filename,wait_flag,your_save_path)

summary_table = readtable('Interneuron_properties_ratio.xlsx');
filenames = cellfun(@(x,y) [x,y],summary_table.Date,summary_table.HFS,'un',0);

file_idx = find(cell2mat(cellfun(@(x) contains(x,filename),filenames,'un',0)));
param.pyr_HS = find(cell2mat(cellfun(@(x) strcmp(x,'Pyr'),{summary_table.HS1{file_idx}, summary_table.HS2{file_idx}},'un',0)));
param.filename = filename;

[d,si,~] = abfload([filename,'.abf']);

if ~isnan(summary_table.last_sweep(file_idx))
    if summary_table.last_sweep(file_idx) > 0
        d_crop = d(:,:,1:summary_table.last_sweep(file_idx));
    elseif summary_table.last_sweep(file_idx) < 0
        d_crop = d(:,:,-summary_table.last_sweep(file_idx):end);
    end
else
    d_crop = d;
end

if size(d,2) > 3
    laser_dim = 5;
    if contains(summary_table.HS1{file_idx}, 'Pyr') || isempty(summary_table.HS1{file_idx})
        int_dim = 3;
        pyr_dim = 1;
    else
        int_dim = 1;
        pyr_dim = 3;
    end
else
    laser_dim = 3;
    if contains(summary_table.HS1{file_idx}, 'Pyr')|| isempty(summary_table.HS1{file_idx})
        int_dim = 2;
        pyr_dim = 1;
    else
        int_dim = 1;
        pyr_dim = 2;
    end    
end

threshold = 0.3;
baseline = 0.05;
signal = 0.02;
sampleFreq = 1/(si*10E-7);
param.laser_V = max(max(squeeze(d_crop(:,laser_dim,:))));

norm_laser = squeeze(d_crop(:,laser_dim,:))/floor(param.laser_V);

laser_on = norm_laser>threshold;

stimOn_idx = sum(cumprod(~laser_on));

nTrials = size(stimOn_idx,2);
for trial_i = 1:nTrials
    baseline_I = mean(d_crop((stimOn_idx(trial_i)-baseline*sampleFreq):stimOn_idx(trial_i),[int_dim pyr_dim],:));
    EPSC.trace(:,:,trial_i) = d_crop(stimOn_idx(trial_i):(stimOn_idx(trial_i)+signal*sampleFreq),[int_dim pyr_dim],trial_i)-baseline_I(:,:,trial_i); 
end
EPSC.meanTrace = mean(EPSC.trace,3);

for HS_i = [1,2] 
    [EPSC.minCurrent(HS_i),min_idx] = min(smooth(EPSC.meanTrace(40:100,HS_i)));
    min_idx = min_idx + 40;
    EPSC.timeToPeak(HS_i) = min_idx;
    EPSC.meanTraceLatency(HS_i) = 30+find(smooth(EPSC.meanTrace(30:end,HS_i))<EPSC.minCurrent(HS_i)*0.2,1,'first');
    total_time = 0.02-(min_idx/sampleFreq);
    [EPSC.decay_tau(HS_i),~,EPSC.tau_err(HS_i)] = fit_exp1(EPSC.meanTrace(min_idx:end,HS_i),total_time);
    temp_min = find(smooth(EPSC.meanTrace(:,HS_i))<EPSC.minCurrent(HS_i)*0.8,1,'first');
    EPSC.meanTraceRise(HS_i) = (temp_min - EPSC.meanTraceLatency(HS_i))/sampleFreq;
    
    % jitter calculation
    for trial_i = 1:nTrials
        EPSC.tMin(HS_i,trial_i) = min(EPSC.trace(1:min_idx,HS_i,trial_i));
        if abs(EPSC.tMin(HS_i,trial_i)) > 3*group_STD
            t1 = find(EPSC.trace(:,HS_i,trial_i)<0.2*EPSC.tMin(HS_i,trial_i),1,'first');
            t2 = find(EPSC.trace(:,HS_i,trial_i)<0.8*EPSC.tMin(HS_i,trial_i),1,'first');
            EPSC.latency(HS_i,trial_i) = t1/sampleFreq;             
        else
            EPSC.latency(HS_i,trial_i) = NaN;
        end
    end    
    EPSC.jitter(HS_i) = nanstd(EPSC.latency(HS_i,:));
end

figure;subplot(1,2,1); hold on;
plot(squeeze(EPSC.trace(:,1,1:5)),'Color',[0.9 0.9 0.9]);
plot(smooth(EPSC.meanTrace(:,1)))
plot(diff(smooth(EPSC.meanTrace(:,1))))
vline([40 120]);
title(num2str(EPSC.jitter(1)*1000));
subplot(1,2,2); hold on;
plot(squeeze(EPSC.trace(:,2,1:5)),'Color',[0.9 0.9 0.9]);
plot(smooth(EPSC.meanTrace(:,2)));
plot(diff(smooth(EPSC.meanTrace(:,2))));
vline([40 120]);
title(num2str(EPSC.jitter(2)*1000));

EPSC.Int_Pyr_ratio = mean(EPSC.minCurrent(1)./EPSC.minCurrent(2));
EPSC.normTrace = EPSC.meanTrace./EPSC.minCurrent;

flag = input('preview EPSCs - keep = enter; exclude = 0; number indicates HS you want to change');
%%
save([your_save_path,filename],'EPSC','param');
try
    xlswrite('Interneuron_properties_ratio.xlsx',...
    {[your_save_path,filename]},1,...
    ['L' num2str(file_idx + 1)]);
    xlswrite('Interneuron_properties_ratio.xlsx',...
    {'X'},1,['K' num2str(file_idx + 1)]);
catch
    wait = input('close excel file');
    xlswrite('Interneuron_properties_ratio.xlsx',...
    {[your_save_path,filename]},1,...
    ['L' num2str(file_idx + 1)]);
    xlswrite('Interneuron_properties_ratio.xlsx',...
    {'X'},1,['K' num2str(file_idx + 1)]);
end

if isempty(flag)
    xlswrite('Interneuron_properties_ratio.xlsx',...
    'X',1,['M' num2str(file_idx + 1)]);
else
    xlswrite('Interneuron_properties_ratio.xlsx',...
    {'N/A'},1,['M' num2str(file_idx + 1)]);
end
end
