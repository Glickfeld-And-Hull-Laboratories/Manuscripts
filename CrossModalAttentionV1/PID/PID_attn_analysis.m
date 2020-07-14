cd('C:\Users\ashley\Documents\Repositories\Manuscripts\CrossModalAttentionV1\PID')

clear all; close all

jb_dir = ['Z:\home\ashley\Manuscripts\Attention V1\Matlab Figs\JB Dataset'];
dataGrabDate = '09-Jul-2020';

load(fullfile(jb_dir,dataGrabDate,'attnData_neuronResp_trialOutcome.mat'))

%% preprocess
npcs=15;  %maxumum number of principle components to save
nc=60;    %number of clusters of patterns (kmeans)

tic
D = preprocess4PID(data,npcs,nc);
toc

%% PID 
binaryS = false;
[PID,MI,ENT] = calcPID(D,npcs,binaryS);

%% plot
 plotinfo