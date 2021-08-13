%% generate spikes structure
% this code analyzes the DC steps data to pass to grouped_HVA_interneurons
% measures Rin, tau, AP adapt, Ih, FI, threshold, spike width
% used in Li 2020 for comparing intrinsic properties across HVAs 

function [spikes,param] = JL2021_intrinsic_props(filename,wait_flag,your_read_path,your_save_path)

% store filename in param struct
param.filename = filename;

% load data
summary_table = readtable('Interneuron_properties_combined.xlsx');

% get abf filenames, load data
filenames = cellfun(@(x,y) [x,y],summary_table.Date,summary_table.DC_Steps,'un',0);
file_idx = find(cell2mat(cellfun(@(x) contains(x,filename),filenames,'un',0)));
[d,si,~] = abfload([filename '.abf']);

% establish analysis params
spikeWindow_S = 0.25;
sampleFreq = 1/(si*10^-6);
delay_S = 0.002;
Vm_analysis_S = 0.1;

% identify which HS is pyr, which HS is int
param.pyr_HS = find(cell2mat(cellfun(@(x) strcmp(x,'Pyr'),{summary_table.HS1{file_idx}, summary_table.HS2{file_idx}},'un',0)));
param.int_HS = find(cell2mat(cellfun(@(x) ~strcmp(x,'Pyr')&&~isempty(x),{summary_table.HS1{file_idx}, summary_table.HS2{file_idx}},'un',0)));

% check if data was recorded for this HS (could be bad cell, cell died
% halfway)
has_data = cell2mat(cellfun(@(x) ~isempty(x),{summary_table.HS1{file_idx} , summary_table.HS2{file_idx}},'un',0));


%% find current injection start and current injection values, store in param struct

Iclamp_HS = find(has_data,1,'first')*2; % if HS1 has data, use line 2 in data matrix, if HS2 has data, use line 4 in data matrix

% get indices 
[~,trace_i] = max(mean(d(:,Iclamp_HS,:)));
Icmd_filt = diff(squeeze(d(:,Iclamp_HS,trace_i)));
[~,on_idx] = max(Icmd_filt);
[~,off_idx] = min(Icmd_filt);

FR_idx = floor(on_idx+(0.2*sampleFreq):on_idx+(0.6*sampleFreq)); % 400 ms 
early_idx = floor(on_idx:on_idx+Vm_analysis_S*sampleFreq);
ss_idx = floor(off_idx-Vm_analysis_S*sampleFreq:off_idx); % 10 ms window for baseline 
rb_idx = floor(off_idx:(off_idx+2*Vm_analysis_S*sampleFreq));

% baseline Vm for future IV analysis
base_IV = squeeze(mean(d(floor(on_idx-(Vm_analysis_S)*sampleFreq:on_idx),:,:)));
baseSub_IV = d-shiftdim(repmat(base_IV,[1,1,size(d,1)]),2);

% current injection values
temp_IV_ss = squeeze(mean(baseSub_IV(ss_idx,:,:)));
param.I_inj = round(temp_IV_ss(Iclamp_HS,:),-1);
param.I_vals = unique(param.I_inj);
param.trial_i = arrayfun(@(x) find(param.I_inj == x),param.I_vals,'un',0);
%% analyze Vm to get spike width, AP adapt, Rin, tau, FI and store in spikes struct

for HS_i = 1:2
    if has_data(HS_i)
        data_dim = HS_i*2 - 1;
        spikes(HS_i).base_IV = base_IV(data_dim,:);
        for sweep_i = 1:size(d,3)
            % late window hyperpolarization analysis
            t_peakVm(sweep_i) = sign(param.I_inj(sweep_i))*max(abs(baseSub_IV(early_idx,data_dim,sweep_i)));
            t_Ih_Vm(sweep_i) = t_peakVm(sweep_i) - min(baseSub_IV(early_idx(end):early_idx(end)+2000,data_dim,sweep_i));
            t_ssVm(sweep_i) = mean(temp_IV_ss(data_dim,sweep_i));
            % get spike times
            spike_idx = on_idx+find(diff(d(on_idx:off_idx,data_dim,sweep_i)>0)==1);
            spikes(HS_i).tSpike_times{sweep_i} = spike_idx/sampleFreq;
            
            if ~isempty(spikes(HS_i).tSpike_times{sweep_i})
                has_spikes(sweep_i) = 1;
                spike_dvdt = diff(d(early_idx(750):end,data_dim,sweep_i));
                threshold_Vm = d(early_idx(750)+find(spike_dvdt > 0.05 * max(spike_dvdt),1,'first'),data_dim,sweep_i);
                spikes(HS_i).tFirst_spike_threshold_Vm(sweep_i) = threshold_Vm;

                [spike_max,spike_idxMax] = max(d(spike_idx(1):spike_idx(1)+(0.001 * sampleFreq),data_dim,sweep_i));
                spike_amp = spike_max-threshold_Vm;
                spike_width_LHS = spike_idx(1)+find(d(spike_idx(1):end,data_dim,sweep_i)>threshold_Vm+(0.5 * spike_amp),1,'first');
                spike_width = find(d(spike_width_LHS:end,data_dim,sweep_i)<threshold_Vm+(0.5 * spike_amp),1,'first');
                spikes(HS_i).tFirst_spike_half_width_ms(sweep_i) = 1000*spike_width/sampleFreq;
                spike_trace(sweep_i,:) = d((spike_idxMax+spike_idx(1)-(0.001 * sampleFreq):spike_idxMax+spike_idx(1)+(0.003*sampleFreq)),data_dim,sweep_i)-spike_max;
                lateWindowSpikes = numel(spike_idx(spike_idx>off_idx-(sampleFreq*spikeWindow_S)));
                earlyWindowSpikes = numel(spike_idx(spike_idx<on_idx+(sampleFreq*0.1)));
                
                if ((lateWindowSpikes == 0) && (earlyWindowSpikes >= 5)) || ((lateWindowSpikes == 0) && isnan(spikes(HS_i).tFR(sweep_i-1)))
                    spikes(HS_i).tFR(sweep_i) = NaN; % depolarization block instances
                    spikes(HS_i).tSpike_rate_adapt(sweep_i) = NaN;
                else
                    subset_idx = (spike_idx>on_idx).*(spike_idx<(on_idx+(2*spikeWindow_S*sampleFreq)));
                    
                    spikes(HS_i).tFR(sweep_i) = sum(subset_idx)/(2*spikeWindow_S); % steady state FR
                    spike_rate = 1./diff(spike_idx(logical(subset_idx))/sampleFreq);
                    
                    if spikes(HS_i).tFR(sweep_i) > 5 
                        spikes(HS_i).tSpike_rate_adapt(sweep_i) = spike_rate(end)/spike_rate(1);
                    else
                        spikes(HS_i).tSpike_rate_adapt(sweep_i) = NaN;
                    end
                end
                
                t_peakVm(sweep_i) = NaN;
                t_ssVm(sweep_i) = NaN;
            else
                has_spikes(sweep_i) = 0;
                spikes(HS_i).tFirst_spike_threshold_Vm(sweep_i) = NaN;
                spikes(HS_i).tFirst_spike_half_width_ms(sweep_i) = NaN;
                spikes(HS_i).tSpike_rate_adapt(sweep_i) = NaN;
                spikes(HS_i).tFR(sweep_i) = 0;
                spike_trace(sweep_i,:) = NaN(1,floor((0.004*sampleFreq)+1));

                [tau_guess(sweep_i)] = fit_exp1(d(on_idx:on_idx+(sampleFreq*0.2),data_dim,sweep_i),0.4);
                
            end
        end
        
        for i = 1:numel(param.I_vals)
            spikes(HS_i).ss_Vm(i) = nanmean(t_ssVm(param.trial_i{i}));
            spikes(HS_i).rebound_Vm(i) = nanmean(t_reboundVm(param.trial_i{i}));
            spikes(HS_i).early_Vm(i) = nanmean(t_peakVm(param.trial_i{i}));
            spikes(HS_i).FR(i) = nanmean(spikes(HS_i).tFR(param.trial_i{i}));
            
            if all(has_spikes(param.trial_i{i}))
                spikes(HS_i).half_width(i) = mean(spikes(HS_i).tFirst_spike_half_width_ms(param.trial_i{i}));
                spikes(HS_i).threshold_Vm(i) = mean(spikes(HS_i).tFirst_spike_threshold_Vm(param.trial_i{i}));
                spikes(HS_i).AP_adapt(i) = nanmean(spikes(HS_i).tSpike_rate_adapt(param.trial_i{i}));
                spikes(HS_i).firstSpike(i,:) = nanmean(spike_trace(param.trial_i{i},:),1);
            else
                spikes(HS_i).half_width(i) = NaN;
                spikes(HS_i).threshold_Vm(i) = NaN;
                spikes(HS_i).AP_adapt(i) = NaN;
                spikes(HS_i).firstSpike(i,:) = NaN(1,size(spike_trace,2));
            end
            
            if sign(param.I_vals(i))<0 % only on hyperpolarizing steps
                Rin_CH(i) = 1000*(abs(spikes(HS_i).ss_Vm(i)./param.I_vals(i)));
                spikes(HS_i).tau_63(i) = nanmean(tau_guess(param.trial_i{i}));
                spikes(HS_i).Ih_Vm(i) = nanmean(t_Ih_Vm(param.trial_i{i}));
            else
                Rin_CH(i) = NaN;
                spikes(HS_i).tau_63(i) = NaN;
                spikes(HS_i).Ih_Vm(i) = NaN;
            end
            
        end
        spikes(HS_i).Rin_CH = nanmean(Rin_CH(~isinf(Rin_CH)));

        %%
        temp_rev = find(param.I_vals>0,1,'first');
        fit_idx = temp_rev-2:temp_rev+1;

        while any(isnan(spikes(HS_i).ss_Vm(fit_idx)))
            temp_idx = fit_idx -1;
            if any(temp_idx<=0)
                fit_idx = temp_idx(2:end);
            else
                fit_idx = temp_idx;
            end
        end
        fit_out = fit(param.I_vals(fit_idx)',spikes(HS_i).ss_Vm(fit_idx)','poly1');
        spikes(HS_i).Rin = fit_out.p1 * 1000;
        
        fit_out = fit(param.I_vals(fit_idx)',spikes(HS_i).early_Vm(fit_idx)','poly1');
        spikes(HS_i).Rin_early = fit_out.p1 * 1000;
        
    else
        
        field_names = [{'base_IV','tSpike_times','tFirst_spike_threshold_Vm','tFirst_spike_half_width_ms','tSpike_rate_adapt','tFR',...
            'ss_Vm','early_Vm','FR','half_width','threshold_Vm','AP_adapt','firstSpike','Ih_Vm','Rin_CH','tau_63','Rin','Rin_early'}];
        
        spikes(HS_i) = cell2struct(mat2cell(NaN(1,numel(field_names)),1,ones(1,numel(field_names))),field_names,2);
        
    end
end
%%

if any([spikes.tau_63]>70)
    wait_flag = true;
end

if wait_flag 
    figure;        
    for HS_i = 1:2
        if has_data(HS_i)
        subplot(2,2,HS_i);
        timeVector = make_time(d,sampleFreq,1);
        plot(timeVector,squeeze(d(:,HS_i*2 - 1,:)));
        axis tight;
        disp(['HS',num2str(HS_i), ' Rin: ',num2str(spikes(HS_i).Rin)]);
        disp(['HS',num2str(HS_i), ' AP adapt: ',num2str(spikes(HS_i).AP_adapt)]);
        
        subplot(2,2,HS_i+2);
        plot(param.I_vals,spikes(HS_i).FR);
        disp(['HS',num2str(HS_i), ' FR: ',num2str(spikes(HS_i).FR)]);
        axis tight;
        
        end
    end

    flag = input('preview spikes - keep = enter; exclude = 0');
else
    flag = [];
end
save([your_save_path,param.filename],'spikes','param');

try
    xlswrite([your_read_path,'Interneuron_properties_combined.xlsx'],...
    {[your_save_path,filename]},1,...
    ['M' num2str(file_idx)]);
catch
    wait = input('close excel file');
    xlswrite([your_read_path,'Interneuron_properties_combined.xlsx'],...
    {[your_save_path,filename]},1,...
    ['M' num2str(file_idx)]);
end
xlswrite([your_read_path,'Interneuron_properties_combined.xlsx'],...
'X',1,['L' num2str(file_idx)]);
if isempty(flag)
    xlswrite([your_read_path,'Interneuron_properties_combined.xlsx'],...
    'X',1,['N' num2str(file_idx)]);
    xlswrite([your_read_path,'Interneuron_properties_combined.xlsx'],...
    'X',1,['O' num2str(file_idx)]);
elseif flag == 0
    xlswrite([your_read_path,'Interneuron_properties_combined.xlsx'],...
    {'N/A'},1,['N' num2str(file_idx)]);
    xlswrite([your_read_path,'Interneuron_properties_combined.xlsx'],...
    {'N/A'},1,['O' num2str(file_idx)]);
elseif flag == 1
    xlswrite([your_read_path,'Interneuron_properties_combined.xlsx'],...
    {'N/A'},1,['N' num2str(file_idx)]);

    xlswrite([your_read_path,'Interneuron_properties_combined.xlsx'],...
    'X',1,['O' num2str(file_idx)]);

elseif flag == 2
    xlswrite([your_read_path,'Interneuron_properties_combined.xlsx'],...
    {'N/A'},1,['O' num2str(file_idx)]);

    xlswrite([your_read_path,'Interneuron_properties_combined.xlsx'],...
    'X',1,['N' num2str(file_idx)]);
end

close all
end
