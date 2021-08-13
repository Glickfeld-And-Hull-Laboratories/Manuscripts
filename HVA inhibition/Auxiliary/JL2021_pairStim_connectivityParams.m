function JL2021_pairStim_connectivityParams(index,u_input,your_read_path,your_save_path)

%% load files
summary_table = readtable([your_read_path,'HVA_IN_pyr_pair_connectivity.xlsx'],'Sheet',u_input);
param.rec_date = summary_table.Date(index);
temp_date = strsplit(summary_table.Date{index},'_');

% addpath(genpath(['Z:\home\jen\Raw data\HVA recordings\',folder_name{u_input},'\',temp_date{1}]));
[d1,si,~] = abfload([temp_date{1},summary_table{index,7}{1},'.abf']);
[d2,~,~] = abfload([temp_date{1},summary_table{index,8}{1},'.abf']);

param.pyrHS = abs(strcmp(summary_table.HS1(index),'Pyr')-2);
param.inHS = abs(param.pyrHS-3);
param.threshold = 0.3;
param.nStims = 10;
param.preStim = 0.01;
param.postStim = 0.04;
param.sampleFreq = si/(10E-4);

data_order = [3,1];
if size(d1,2) == 4
    stim_order = [2,4];
elseif size(d1,2) == 3
    stim_order = [2,3];
end

if param.pyrHS == 1
    pyr_in_data = squeeze(d1(:,data_order(1),:));
    in_pyr_data = squeeze(d2(:,data_order(2),:));
    pyr_in_stim = squeeze(d1(:,stim_order(1),:));
    in_pyr_stim = squeeze(d2(:,stim_order(2),:));
elseif param.pyrHS == 2
    pyr_in_data = squeeze(d2(:,data_order(2),:));
    in_pyr_data = squeeze(d1(:,data_order(1),:));
    pyr_in_stim = squeeze(d2(:,stim_order(2),:));
    in_pyr_stim = squeeze(d1(:,stim_order(1),:));
end


%% get average traces

% Pyr-->IN

stimOn = find(diff(pyr_in_stim>(param.threshold*max(pyr_in_stim))>0));
stimOn = stimOn(1:2:end);

% if no stim On found
if isempty(stimOn)
    % check other file
    stimOn = find(diff(in_pyr_stim>(param.threshold*max(in_pyr_stim))>0));
    stimOn = stimOn(1:2:end);
end

crop_data = arrayfun(@(x) pyr_in_data(x-(param.preStim*param.sampleFreq):x+(param.postStim*param.sampleFreq))-...
    mean(pyr_in_data(x-(2*param.preStim*param.sampleFreq):x-(param.preStim*param.sampleFreq))),stimOn,'un',0);

for stim = 1:param.nStims 
    data.IN_EPSC_all{stim} = cell2mat(crop_data(stim:10:end));
end

data.IN_EPSC_mean = cellfun(@(x) mean(x),data.IN_EPSC_all,'un',0);
data.IN_EPSC_amp = cell2mat(cellfun(@(x) min(smooth(x(param.preStim*param.sampleFreq:end))),data.IN_EPSC_mean,'un',0));

% IN --> Pyr

stimOn = find(diff(in_pyr_stim>(param.threshold*max(in_pyr_stim))>0));
stimOn = stimOn(1:2:end);

if isempty(stimOn)
    stimOn = find(diff(pyr_in_stim>(param.threshold*max(pyr_in_stim))>0));
    stimOn = stimOn(1:2:end);
end

crop_data = arrayfun(@(x) in_pyr_data(x-(param.preStim*param.sampleFreq):x+(param.postStim*param.sampleFreq))-...
    mean(in_pyr_data(x-(2*param.preStim*param.sampleFreq):x-(param.preStim*param.sampleFreq))),stimOn,'un',0);

for stim = 1:param.nStims 
    data.pyr_IPSC_all{stim} = cell2mat(crop_data(stim:param.nStims:end));
end

data.pyr_IPSC_mean = cellfun(@(x) mean(x),data.pyr_IPSC_all,'un',0);
data.pyr_IPSC_amp = cell2mat(cellfun(@(x) max(smooth(x(param.preStim*param.sampleFreq:end))),data.pyr_IPSC_mean,'un',0));


savePathName = [your_save_path,param.rec_date{1}];

if ~exist(savePathName,'dir')
    mkdir(savePathName);
    addpath(savePathName)
    tempPath = strcat(savePathName,'\pathdef.m');
    savepath tempPath
end


try
    xlswrite([your_read_path,'HVA_IN_pyr_pair_connectivity.xlsx'],{[savePathName,'\',param.rec_date{1},'_pairConnectivity']},u_input,['K' num2str(index+1)]);
catch
    input('close the excel file. press enter when done.');
    xlswrite([your_read_path,'HVA_IN_pyr_pair_connectivity.xlsx'],{[savePathName,'\',param.rec_date{1},'_pairConnectivity']},u_input,['K' num2str(index+1)]);
end
   
save([savePathName,'\',param.rec_date{1},'_pairConnectivity'],'data','param');
