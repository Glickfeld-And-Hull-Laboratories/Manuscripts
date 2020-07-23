cd('C:\Users\ashley\Documents\Repositories\Manuscripts\CrossModalAttentionV1\PID')

clear all; close all

%% load previous?
runTag = 'boots_exptOr15npcs_categorical';
doLoadPrevious = false;
doBootstrapPID = true;
doMultiPC = false;
binaryS = false;
nBoot = 100;
% prevTimestampID = '200716_1652'; %not bootstrapped, expt npcs, binary
% prevTimestampID = '200721_1946'; %not bootstrapped, expt npcs, categorical
% prevTimestampID = '200716_1841'; %not bootstrapped, 15 npcs
% prevTimestampID = '200720_2003'; %'200716_1436';% boostrapped
prevTimestampID = '200720_1802'; %not bootstrapped, multiple npcs
%%

jb_dir = ['Z:\home\ashley\Manuscripts\Attention V1\Matlab Figs\JB Dataset'];
dataGrabDate = '09-Jul-2020';

load(fullfile(jb_dir,dataGrabDate,'attnData_neuronResp_trialOutcome.mat'))

fnout = fullfile(jb_dir,'PID Analysis');
% mkdir(fnout)
%% preprocess
if doMultiPC 
    npcs = [5,10,15,20];    %maxumum number of principle components to save
    nc=60;    %number of clusters of patterns (kmeans)
    d = cell(1,length(npcs));
    for ipc = 1:length(npcs)
        tic
        D{ipc} = preprocess4PID(data,npcs(ipc),nc,binaryS);
        t=toc; 
    end
else
    npcs=15;  %maxumum number of principle components to save
    nc=60;    %number of clusters of patterns (kmeans)

    tic
    D = preprocess4PID(data,'expt_plus',nc,binaryS);
    t=toc;    
end
%% PID 
if doLoadPrevious
    load(fullfile(fnout,['PID_results_' prevTimestampID]))
elseif doBootstrapPID
    PIDresults = struct;
    tic
    for iboot = 1:nBoot
        D_boot = cell(1,length(D));
        for iexp = 1:length(D) %grab random trials for each experiment
            nt = [length(D{iexp}.Vis.B),length(D{iexp}.Aud.B)];
            rand_trial_idx = arrayfun(@(x) randsample(x,x,0),nt,'unif',0);
            d = D{iexp};
            for iav = 1:2
                if iav == 1
                    dd = d.Vis;
                else
                    dd = d.Aud;
                end
                dd.Rlabels = dd.Rlabels(rand_trial_idx{iav},:);
                dd.Rraw = dd.Rraw(rand_trial_idx{iav},:);
                dd.B = dd.B(rand_trial_idx{iav},:);
                dd.R = dd.R(rand_trial_idx{iav},:);
                dd.S = dd.S(rand_trial_idx{iav},:);
                dd.NS = dd.NS(rand_trial_idx{iav},:);
                dd.frame = dd.frame(rand_trial_idx{iav},:);
                if iav == 1
                    d.Vis = dd;
                else
                    dd.Aud = dd;
                end                
            end
            D_boot{iexp} = d;
        end
        fprintf('Starting PID bootstrap %s...\n',num2str(iboot))
        [PID,MI,ENT] = calcPID(D_boot,npcs,binaryS);
        PIDresults.PID{iboot} = PID;
        PIDresults.MI{iboot} = MI;
        PIDresults.ENT{iboot} = ENT;
    end
    PIDresults.binaryS = binaryS;
    PIDresults.D = D;
    toc
    timestampID = [datestr(now,'yymmdd') '_' datestr(now,'HHMM')];
    save(fullfile(fnout,['PID_results_' timestampID]),'PIDresults')
elseif doMultiPC
    PIDresults = struct;
    for ipc = 1:length(npcs)
        [PID,MI,ENT] = calcPID(D{ipc},npcs(ipc),binaryS);
        PIDresults.PID{ipc} = PID;
    end
    PIDresults.binaryS = binaryS;
    PIDresults.D = D;
    timestampID = [datestr(now,'yymmdd') '_' datestr(now,'HHMM')];
    save(fullfile(fnout,['PID_results_' timestampID]),'PIDresults')
else
    tic
    [PID,MI,ENT] = calcPID(D,npcs,binaryS);
    toc
    PIDresults.PID = PID;
    PIDresults.MI = MI;
    PIDresults.ENT = ENT;
    PIDresults.binaryS = binaryS;
    PIDresults.D = D;
    timestampID = [datestr(now,'yymmdd') '_' datestr(now,'HHMM')];
    save(fullfile(fnout,['PID_results_' timestampID]),'PIDresults')
end
%% plot
plotinfo_aw