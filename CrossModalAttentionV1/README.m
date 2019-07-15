%% Below are the major analyses performed for the manuscirpt "Separable codes for read-out of mouse primary visual cortex across attentional states"
% Several functions used in these analysis are not listed, but are part of
% the ImagingCode-Glickfeld-Hull or BehaviorCode-Glickfeld-Hull
% repositories in the "Hull and Glickfeld Laboratories" organiation's
% GitHub page.

%% Behavior data analysis
% eaMsBxSummary_attnV1ms - process each mouse's behavior sessions
% bxSummary_FSAV_attnV1ms - analyze data and make figs for manuscript
% bxParams_FSAV_attnV1ms - parameters, such as behavior cutoffs, used for analysis
%% Pupil data analysis
% getPupilTC - extract size and position from eye movies
% createFSAVEyeStruct - organize data for anticipation and target aligned analyses
% exampleEyeImage_FSAVattentionV1_attnV1ms - figure for example eye image
% plotFSAV_eye_attnV1ms - analyze data and make figures 
%% Imaging data analysis
% FSAV_attentionV1 - behavior dataset info
% FSAV_V1_100ms_naive - naive dataset info
% compareRegImg - get registration image
% reRegDatasets - register data set and select ROIs
% getOriTuningFSAV - get orienation tuning curves and responses 
% createFSAVDataStruct - organize data for anticipation and target aligned analyses - enter in behavior or naive dataset name
% createFSAVDataStruct2 - organize data for decoding models - enter in behavior or naive dataset name
% FSAV_imgAnalysis_anticipation_attnV1ms - analyze data and make figs for manuscript - enter in behavior or naive dataset name