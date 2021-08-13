%% Supplementary figure 5

%% Panels A-B
control_name = {'20d28008';'20d28024';'20d28037'};
muscimol_name = {'20d28009';'20d28025';'20d28038'};
ttx_name = {'20d28012';'20d28026';'20d28039'};

for exp_i = 1:numel(control_name)
control_dat = abfload([control_name{exp_i} '.abf']);
muscimol_dat = abfload([muscimol_name{exp_i} '.abf']);
ttx_dat = abfload([ttx_name{exp_i} '.abf']);

stimOn_idx = repelem(8313,size(control_dat,3));
LFP.control = cleanData(control_dat,1,stimOn_idx,'afterStimS',0.01);
LFP.muscimol = cleanData(muscimol_dat,1,stimOn_idx,'afterStimS',0.01);
LFP.TTX = cleanData(ttx_dat,1,stimOn_idx,'afterStimS',0.01);

figure; hold on;
plot(timeVector,mean(LFP.control.pts,2)-mean(LFP.TTX.pts,2));
plot(timeVector,mean(LFP.muscimol.pts,2)-mean(LFP.TTX.pts,2));
fix_axes(gcf,20,'time','mV');
end

%% Panel C
threshold = 0.3;
baseline = 0.05;
signal = 0.02;
sampleFreq = 20000;

% load CGP
temp_dat = abfload('20d28015.abf');
temp_subDat = squeeze(temp_dat(:,2,2:2:end));
norm_laser = temp_subDat/5;
laser_on = norm_laser>threshold;
stimOn_idx = sum(cumprod(~laser_on));

temp_dat = temp_dat(:,1,2:2:end);
CGP_dat = cleanData(temp_dat,1,stimOn_idx,'afterStimS',0.015);

% load + 12.5 musc
temp_dat = abfload('20d28016.abf');
temp_dat = temp_dat(:,1,2:2:end);
muscLow_dat = cleanData(temp_dat,1,stimOn_idx,'afterStimS',0.015);

% load + 25 musc
temp_dat = abfload('20d28019.abf');
temp_dat = temp_dat(:,1,2:2:end);
muscMed_dat = cleanData(temp_dat,1,stimOn_idx,'afterStimS',0.015);

% load +50 musc
temp_dat = abfload('20d28021.abf');
temp_dat = temp_dat(:,1,2:2:end);
muscHi_dat = cleanData(temp_dat,1,stimOn_idx,'afterStimS',0.015);

timeVector = make_time(CGP_dat.pts,20,1);
figure; subplot(1,2,1);
hold on;
plot(timeVector,smooth(mean(CGP_dat.pts,2)));
plot(timeVector,mean(muscLow_dat.pts,2));
plot(timeVector,mean(muscMed_dat.pts,2));
plot(timeVector,smooth(mean(muscHi_dat.pts,2)));
legend({'+CGP';'+10 uM musc';'+25 uM musc';'+50 um musc'},'Location','northoutside');
axis square;

title('raw LFP - 5 V');
fix_axes(gcf,10,'time (ms)','mV');
