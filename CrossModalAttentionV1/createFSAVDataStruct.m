function mouse = createFSAVDataStruct(datasetStr,cellsOrDendrites,varargin)
%cellsOrDendrites: 1 == cells; 2 == dendrites; 3 == axons
if nargin > 2
    doRawData = varargin{1};
else
    doRawData = 0;
end
    % set analysis windows
    pre_event_time = 1000; %ms
    post_event_time = 4500; %ms
    pre_win_time = [-30 70];
    trans_win_time = [150 250]; %this is actually ~166-266 ms at 30 Hz
    minTrialLengthMs = 2500; % time in ms of trial rather than hard-coding cycle number for analysis
    
    
    rc = behavConstsAV;

    eval(datasetStr)

    dataGroup = datasetStr;
    
    %create list of mice for file names
    nMice = length(unique({expt.SubNum}));
    str = unique({expt.SubNum});
    values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
    mouse_str = ['i' strjoin(str,'_i')];
    %initialize structure
    s = zeros(1,nMice);
    mouse = struct;
    Dirs_all = [];
    
    setFigParams4Print('portrait')
    nexp = size(expt,2);
    n = ceil(sqrt(nexp+1));
    if (n^2)-n > nexp+1
        n2 = n-1;
    else
        n2 = n;
    end
    motionHist = figure;
    
    %collect data from each experiment
    for iexp = 1:size(expt,2)
        disp([num2str(expt(iexp).date) ' i' num2str(expt(iexp).SubNum)])

        nrun = size(expt(iexp).runs,1);
        dir_run = expt(iexp).dirtuning;
        
        %keep track of which mouse
        imouse = find(values == str2num(expt(iexp).SubNum));
        s(:,imouse) = s(:,imouse)+1;
        mouse(imouse).expt(s(:,imouse)).date = expt(iexp).date;
        mouse(imouse).expt(s(:,imouse)).mouse_name = expt(iexp).SubNum;
        
        %translate time windows into frames 
        pre_event_frames = ceil(pre_event_time*(expt(iexp).frame_rate/1000));
        post_event_frames = ceil(post_event_time*(expt(iexp).frame_rate/1000));
        pre_win_frames = pre_event_frames+round(pre_win_time.*(expt(iexp).frame_rate/1000));
        pre_win = pre_win_frames(1):pre_win_frames(2);
        trans_win_frames = pre_event_frames+round(trans_win_time.*(expt(iexp).frame_rate/1000));
        trans_win = trans_win_frames(1):trans_win_frames(2);
        minTrialLengthFrames = ceil(minTrialLengthMs*(expt(iexp).frame_rate/1000));
                
        %create string for saving mult runs
        runstr = expt(iexp).runs(1,:);
        if nrun>1
            for irun = 2:nrun
                runstr = [runstr '-' expt(iexp).runs(irun,:)];
            end
        end
        fnout = fullfile(rc.caOutputDir, expt(iexp).mouse, expt(iexp).folder, expt(iexp).date, [expt(iexp).date '_' expt(iexp).mouse '_' runstr '_']);
        
        % load and combine mworks data and timecourses
        input = [];
        for irun = 1:nrun
            time = expt(iexp).time_mat(irun,:);
            fn_mworks = [rc.pathStr...
                '\data-i' expt(iexp).SubNum '-' expt(iexp).date '-' time '.mat'];
            if irun == 1
                input = mwLoadData(fn_mworks, [], []);
            else
                try
                    input = [input mwLoadData(fn_mworks, [], [])];
                catch
                    input2 = mwLoadData(fn_mworks, [], []);
                    inpNames1 = fieldnames(input);
                    inpNames2 = fieldnames(input2);
                    inpLong = gt(length(inpNames1),length(inpNames2));
                    if inpLong == 1
                        inpPlusInd = ismember(inpNames1,inpNames2);
                        inpPlus = inpNames1(~inpPlusInd);
                        for i = 1:length(inpPlus)
                            input2.(genvarname(inpPlus{i})) = ...
                                cell(1,input2.trialSinceReset);
                        end
                    else
                        inpPlusInd = ismember(inpNames2,inpNames1);
                        inpPlus = inpNames2(~inpPlusInd);
                        for i = 1:length(inpPlus)
                            input.(char(genvarname(inpPlus(i)))) = cell(1,80);
                        end
                    end
                    input = [input input2];
                end
            end
        end
        input = concatenateDataBlocks(input);  
        if isnan(expt(iexp).trial_range)
            tr = 1:length(input.trialOutcomeCell);
        else
            tr = expt(iexp).trial_range;
        end     
        
        %convert important fields to matrices
        run_trials = input.trialsSinceReset;
        cLeverDown = celleqel2mat_padded(input.cLeverDown);
        cFirstStim = celleqel2mat_padded(input.cFirstStim);
        cLeverUp = celleqel2mat_padded(input.cLeverUp);
        cTargetOn = celleqel2mat_padded(input.cTargetOn);
        cStimOn = celleqel2mat_padded(input.cStimOn);
        cItiStart = celleqel2mat_padded(input.cItiStart);
        reactTimes = celleqel2mat_padded(input.reactTimesMs);
        reactTimes = reactTimes(tr);
        if expt(iexp).catch
        cCatchOn = celleqel2mat_padded(input.cCatchOn);
        isFA = celleqel2mat_padded(input.tFalseAlarm);
        tCatchGratingDirectionDeg = chop(celleqel2mat_padded(input.tCatchGratingDirectionDeg),4);
        
        isFA = isFA(tr);
        tCatchGratingDirectionDeg(tr);
        
        cDirs = unique(tCatchGratingDirectionDeg(~isnan(tCatchGratingDirectionDeg)));        
        end
        tooFastTime = input.nFramesTooFast;
        maxReactTime = input.nFramesReact;
        


        
        %account for accumulation of frames across multiple runs 
        dataTC = [];
        if strcmp(datasetStr,'awFSAVdatasets_V1axonsAL') | strcmp(datasetStr,'awFSAVdatasets_V1axonsPM')
            fnTC = fullfile(rc.ashleyAnalysis,...
                expt(iexp).mouse,expt(iexp).folder, expt(iexp).date);
        else
            fnTC = fullfile(rc.ashleyAnalysis,...
                expt(iexp).mouse,expt(iexp).folder, expt(iexp).date,'data processing');
        end

        if cellsOrDendrites == 1
                load(fullfile(fnTC,'timecourses_bx_cells.mat'))
            if doRawData == 1
                dataTC = data_bx_tc;
            else
                try
                    dataTC = data_bx_tc_subnp;
                catch
                    dataTC = dataTC_npSub;
                end
            end
        elseif cellsOrDendrites == 2
            load(fullfile(fnTC,'timecourses_bx_dendrites.mat'))
            if doRawData == 1
                dataTC = data_bx_den_tc;
            else
                dataTC = data_bx_den_tc_subnp;
            end
        elseif cellsOrDendrites == 3
            load(fullfile(fnTC,'timecourses.mat'))
            dataTC = data_tc_subnp;
        end
                
        offset = 0;
        for irun = 1:nrun
            ImgFolder = expt(iexp).runs(irun,:);
            nfr_run = nFramesSbxDataset(expt(iexp).mouse,expt(iexp).date,ImgFolder);
            offset = offset+nfr_run;
            if irun < nrun
                startTrial = sum(run_trials(1, 1:irun),2)+1;
                endTrial = sum(run_trials(1,1:irun+1),2);
                cLeverDown(1,startTrial:endTrial) = cLeverDown(1,startTrial:endTrial)+offset;
                cFirstStim(1,startTrial:endTrial) = cFirstStim(1,startTrial:endTrial)+offset;
                cLeverUp(1,startTrial:endTrial) = cLeverUp(1,startTrial:endTrial)+offset;
                cTargetOn(1,startTrial:endTrial) = cTargetOn(1,startTrial:endTrial)+offset;
                cStimOn(1,startTrial:endTrial) = cStimOn(1,startTrial:endTrial)+offset;
                cItiStart(1,startTrial:endTrial) = cItiStart(1,startTrial:endTrial)+offset;
                if expt(iexp).catch
                    cCatchOn(1,startTrial:endTrial) = cCatchOn(1,startTrial:endTrial)+offset;
                end
            end
        end
        
        cLeverDown = cLeverDown(tr);
        cFirstStim = cFirstStim(tr);
        cLeverUp = cLeverUp(tr);
        cTargetOn = cTargetOn(tr);
        cStimOn = cStimOn(tr);
        cItiStart = cItiStart(tr);
        if expt(iexp).catch
            cCatchOn = cCatchOn(tr);
        end

        ntrials = length(tr);%length(input.trialOutcomeCell);
        
        tCyclesOn = cell2mat(input.tCyclesOn);
        minCyclesOn = input.minCyclesOn;
        V_ind = find(cell2mat(input.tBlock2TrialNumber(tr)) == 0);
        AV_ind = find(cell2mat(input.tBlock2TrialNumber(tr)) == 1);
        
        tCyclesOn = tCyclesOn(tr);
        
        %previous trial info
        trType = double(cell2mat(input.tBlock2TrialNumber));
        
        trType = trType(tr);
        
        trType_shift = [NaN trType];
        prevTrType = trType_shift(1:length(trType));
        trType_shift = [NaN NaN trType];
        prev2TrType = trType_shift(1:length(trType));
        
        
        trOut = input.trialOutcomeCell;
        trOut = trOut(tr);
        
        trOut_shift = [{NaN} trOut];
        prevTrOut = trOut_shift(1:length(trOut));
        trOut_shift = [{NaN} {NaN} trOut];
        prev2TrOut = trOut_shift(1:length(trOut));
        
        %stim timing
        if iscell(input.nFramesOn)
            cycTime = unique(cell2mat(input.nFramesOn))+unique(cell2mat(input.nFramesOff));
        else
            cycTime = input.nFramesOn+input.nFramesOff;
        end
        minTime = minCyclesOn*cycTime;
        frameratems = expt(iexp).frame_rate/1000;
        cycTimeMs = cycTime/frameratems;
        minCyclesAnt = floor(minTrialLengthFrames/cycTime);

        % for older datasets that do not have catch trial outcome in input
        %disp(sum(cell2mat(input.tShortCatchTrial)))
        if expt(iexp).catch
            catchIndex = cell2mat(input.tShortCatchTrial);
            catchIndex = find(catchIndex(tr));
            catchReactFrames = cLeverUp-cCatchOn;
            catchTrialReactTimeMs = NaN(1,ntrials);
            catchTrialReactTimeMs(catchIndex) = catchReactFrames(catchIndex)/frameratems;
            catchCyclesOn = cell2mat_padded(input.catchCyclesOn);
            catchCyclesOn = catchCyclesOn(tr);
            
            if ~any(strcmp(fieldnames(input),'catchTrialOutcomeCell'))
                catchTrialOutcome = cell(1,ntrials);
                for i = 1:sum(cell2mat(input.tShortCatchTrial(tr)))
                    if isFA(catchIndex(i)) == 1
                        catchTrialOutcome{catchIndex(i)} = 'FA';
                    elseif cCatchOn(catchIndex(i)) == 0
                        catchTrialOutcome{catchIndex(i)} = 'failure';
                    elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) < tooFastTime
                        catchTrialOutcome{catchIndex(i)} = 'failure';
                    elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) > maxReactTime
                        catchTrialOutcome{catchIndex(i)} = 'CR';
                    end
                end
            else
                if length(input.catchTrialOutcomeCell) < ntrials
                    catchTrialOutcome = cell(1,ntrials);
                    for i = 1:sum(cell2mat(input.tShortCatchTrial(tr)))
                        if isFA(catchIndex(i)) == 1
                            catchTrialOutcome{catchIndex(i)} = 'FA';
                        elseif isnan(cCatchOn(catchIndex(i)))
                            catchTrialOutcome{catchIndex(i)} = 'failure';
                        elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) < tooFastTime
                            catchTrialOutcome{catchIndex(i)} = 'failure';
                        elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) > maxReactTime
                            catchTrialOutcome{catchIndex(i)} = 'CR';
                        end
                    end
                elseif length(input.catchTrialOutcomeCell) < tr(end)
                    trInd = tr(1:find(tr == length(input.catchTrialOutcomeCell)));
                    nPaddedTr = sum(tr > length(input.catchTrialOutcomeCell));
                    catchTrialOutcome = cat(2,...
                        input.catchTrialOutcomeCell(trInd),cell(1,nPaddedTr));
                    
                else
                    catchTrialOutcome = input.catchTrialOutcomeCell(tr);
                end
            end
        end

        tGratingDirectionDeg = chop(celleqel2mat_padded(input.tGratingDirectionDeg),4);
        tGratingDirectionDeg = tGratingDirectionDeg(tr);
        Dirs = unique(tGratingDirectionDeg);
        Dirs_all = unique([Dirs_all Dirs]);
        
        %load direction tuning data
        fnTun = fullfile(rc.ashleyAnalysis,...
            expt(iexp).mouse,expt(iexp).folder, expt(iexp).date, ...
            expt(iexp).dirtuning);
        if cellsOrDendrites == 1
            load(fullfile(fnTun, 'oriTuningAndFits.mat'));
        elseif cellsOrDendrites == 2
            load(fullfile(fnTun, 'oriTuningAndFits_den.mat'));
        elseif cellsOrDendrites == 3
            load(fullfile(fnTun, 'oriTuningAndFits_boutons.mat'));
        end
        
        mouse(imouse).expt(s(:,imouse)).cells.oriTuning.oriResp = avgResponseEaOri;
        mouse(imouse).expt(s(:,imouse)).cells.oriTuning.oriRespSem = semResponseEaOri;
        mouse(imouse).expt(s(:,imouse)).cells.oriTuning.oriFit = vonMisesFitAllCells;
        mouse(imouse).expt(s(:,imouse)).cells.oriTuning.oriFitReliability = ...
            fitReliability;
        
        %sort by trial type
        if expt(iexp).catch
            Ix = find(cell2mat(input.tShortCatchTrial(tr))==0);
            CIx = find(cell2mat(input.tShortCatchTrial(tr))==1);
        else
            Ix = 1:length(input.trialOutcomeCell(tr));
        end
        
        F1Ix = find(tCyclesOn == 1);
        F1oneStim = intersect(V_ind,F1Ix);
        F2oneStim = intersect(AV_ind,F1Ix);
        
        FIx = intersect(Ix, find(strcmp(input.trialOutcomeCell(tr), 'failure')));
        SIx = intersect(intersect(Ix, find(reactTimes>200)), find(strcmp(input.trialOutcomeCell(tr), 'success')));
        SIxAllReact = intersect(Ix, find(strcmp(input.trialOutcomeCell(tr), 'success')));
        MIx = intersect(Ix, find(strcmp(input.trialOutcomeCell(tr), 'ignore')));
        FIxlong = intersect(find(tCyclesOn>3), FIx);
        SIxlong = intersect(find(tCyclesOn>3), SIx);
        MIxlong = intersect(find(tCyclesOn>3), MIx);
        Fb1Ix = intersect(V_ind, FIxlong);
        Fb2Ix = intersect(AV_ind, FIxlong);
        Sb1Ix = intersect(V_ind, SIxlong);
        Sb2Ix = intersect(AV_ind, SIxlong);
        Mb1Ix = intersect(V_ind, MIxlong);
        Mb2Ix = intersect(AV_ind, MIxlong);
        Rb1Ix = intersect(Ix,intersect(V_ind, find(tCyclesOn>3)));
        Rb2Ix = intersect(Ix,intersect(AV_ind, find(tCyclesOn>3)));
        SIxAllReactlong = intersect(find(tCyclesOn>3), SIxAllReact);
        
        %find result of previous trial
        if sum(FIx==1)
            if sum(intersect(FIx,V_ind)==1)
            Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx(2:end)-1,V_ind)));
            Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx(2:end)-1,AV_ind)));
            Fb5Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,V_ind)));
            Fb6Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,AV_ind)));
            else
            Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
            Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
            Fb5Ix = intersect(Fb2Ix,FIx(ismember(FIx(2:end)-1,V_ind)));
            Fb6Ix = intersect(Fb2Ix,FIx(ismember(FIx(2:end)-1,AV_ind)));
            end
            Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
            Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
            Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,V_ind)));
            Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,AV_ind)));
            Sb5Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,V_ind)));
            Sb6Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,AV_ind)));
            Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,V_ind)));
            Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,AV_ind)));
            SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
            SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
            SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
            SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
        elseif sum(SIx==1)            
            if sum(intersect(SIx,V_ind)==1)
            Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx(2:end)-1,V_ind)));
            Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx(2:end)-1,AV_ind)));
            Sb5Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,V_ind)));
            Sb6Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,AV_ind)));
            else
            Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
            Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
            Sb5Ix = intersect(Sb2Ix,SIx(ismember(SIx(2:end)-1,V_ind)));
            Sb6Ix = intersect(Sb2Ix,SIx(ismember(SIx(2:end)-1,AV_ind)));
            end
            Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
            Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
            Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,V_ind)));
            Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,AV_ind)));
            Fb5Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,V_ind)));
            Fb6Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,AV_ind)));
            Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,V_ind)));
            Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,AV_ind)));
            SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
            SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
            SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
            SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
        elseif sum(MIx==1)
            if sum(intersect(MIx,V_ind)==1)
            Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx(2:end)-1,V_ind)));
            Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx(2:end)-1,AV_ind)));
            Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,V_ind)));
            Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,AV_ind)));
            else
            Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,V_ind)));
            Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,AV_ind)));
            Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx(2:end)-1,V_ind)));
            Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx(2:end)-1,AV_ind)));
            end
            Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
            Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
            Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
            Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
            Fb5Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
            Fb6Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
            Sb5Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
            Sb6Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
            SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
            SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
            SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
            SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
        elseif sum(SIxAllReactlong==1)
            if sum(intersect(SIxAllReactlong,V_ind)==1)
            SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong(2:end)-1,V_ind)));
            SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong(2:end)-1,AV_ind)));
            SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
            SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
            else
            SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
            SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
            SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong(2:end)-1,V_ind)));
            SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong(2:end)-1,AV_ind)));                
            end
            Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
            Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
            Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
            Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
            Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,V_ind)));
            Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,AV_ind)));
            Fb5Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,V_ind)));
            Fb6Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,AV_ind)));
            Sb5Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,V_ind)));
            Sb6Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,AV_ind)));
            Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,V_ind)));
            Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,AV_ind)));
        else
            Fb3Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,V_ind)));
            Fb4Ix = intersect(Fb1Ix,FIx(ismember(FIx-1,AV_ind)));
            Sb3Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,V_ind)));
            Sb4Ix = intersect(Sb1Ix,SIx(ismember(SIx-1,AV_ind)));
            Mb3Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,V_ind)));
            Mb4Ix = intersect(Mb1Ix,MIx(ismember(MIx-1,AV_ind)));
            Fb5Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,V_ind)));
            Fb6Ix = intersect(Fb2Ix,FIx(ismember(FIx-1,AV_ind)));
            Sb5Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,V_ind)));
            Sb6Ix = intersect(Sb2Ix,SIx(ismember(SIx-1,AV_ind)));
            Mb5Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,V_ind)));
            Mb6Ix = intersect(Mb2Ix,MIx(ismember(MIx-1,AV_ind)));
            SbAR3Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
            SbAR4Ix = intersect(V_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
            SbAR5Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,V_ind)));
            SbAR6Ix = intersect(AV_ind,SIxAllReactlong(ismember(SIxAllReactlong-1,AV_ind)));
        end
            
        
        SbAR1Ix = intersect(V_ind, SIxAllReactlong);
        SbAR2Ix = intersect(AV_ind, SIxAllReactlong);

        %find direction with maximum hits and misses to match datasets
        Sb1IxMatch = [];
        Mb1IxMatch = [];
        nS = zeros(1,length(Dirs));
        nM = zeros(1,length(Dirs));
        for iDir = 1:length(Dirs)
            dirIx = find(tGratingDirectionDeg==Dirs(iDir));
            Sb1IxDir = intersect(dirIx,Sb1Ix);
            Mb1IxDir = intersect(dirIx,Mb1Ix);
            nS(iDir) = length(Sb1IxDir);
            nM(iDir) = length(Mb1IxDir);
            if nS(iDir)<nM(iDir)
                Sb1IxMatch = [Sb1IxMatch Sb1IxDir];
                Mb1IxMatch = [Mb1IxMatch Mb1IxDir(randperm(nM(iDir),nS(iDir)))];
            elseif nS(iDir)>nM(iDir)
                Mb1IxMatch = [Mb1IxMatch Mb1IxDir];
                Sb1IxMatch = [Sb1IxMatch Sb1IxDir(randperm(nS(iDir),nM(iDir)))];
            else
                Sb1IxMatch = [Sb1IxMatch Sb1IxDir];
                Mb1IxMatch = [Mb1IxMatch Mb1IxDir];
            end
        end
        if isfield(input, 'tSoundTargetAmplitude')
            tSoundTargetAmplitude = celleqel2mat_padded(input.tSoundTargetAmplitude(tr))';
            Amps = unique(chop(tSoundTargetAmplitude,2));
            nS = zeros(1,length(Amps));
            nM = zeros(1,length(Amps));
            Sb2IxMatch = [];
            Mb2IxMatch = [];
            if length(Amps)>2
                for iAmp = 2:length(Amps)
                    ampIx = find(tSoundTargetAmplitude==Amps(iAmp));
                    Sb2IxAmp = intersect(ampIx,Sb2Ix);
                    Mb2IxAmp = intersect(ampIx,Mb2Ix);
                    nS(iAmp) = length(Sb2IxAmp);
                    nM(iAmp) = length(Mb2IxAmp);
                    if (nS(iAmp)+nM(iAmp))>0
                        if nS(iAmp)<nM(iAmp)
                            Sb2IxMatch = [Sb2IxMatch; Sb2IxAmp];
                            Mb2IxMatch = [Mb2IxMatch; Mb2IxAmp(randperm(nM(iAmp),nS(iAmp)))];
                        elseif nS(iAmp)>nM(iAmp)
                            Mb2IxMatch = [Mb2IxMatch; Mb2IxAmp];
                            Sb2IxMatch = [Sb2IxMatch; Sb2IxAmp(randperm(nS(iAmp),nM(iAmp)))];
                        else
                            Sb2IxMatch = [Sb2IxMatch; Sb2IxAmp];
                            Mb2IxMatch = [Mb2IxMatch; Mb2IxAmp];
                        end
                    end
                end
            else
                Sb2IxMatch = Sb2Ix;
                Mb2IxMatch = Mb2Ix;
            end
        else
            Amps = [0 input.soundTargetAmplitude];
            Sb2IxMatch = Sb2Ix;
            Mb2IxMatch = Mb2Ix;
        end
        
    % match number of trials per direction for invalid and valid hits and
    % misses (catch trials only)
    if expt(iexp).catch
        fa1Ix = find(strcmp(catchTrialOutcome,'FA'));
        cr1Ix = find(strcmp(catchTrialOutcome,'CR'));
        fa1Ix = fa1Ix;
        cr1Ix = cr1Ix;
        if sum(fa1Ix==1)
            fa3Ix = intersect(fa1Ix, fa1Ix(ismember(fa1Ix(2:end)-1,V_ind)));
            fa4Ix = intersect(fa1Ix, fa1Ix(ismember(fa1Ix(2:end)-1,AV_ind)));
            cr3Ix = intersect(cr1Ix, cr1Ix(ismember(cr1Ix-1,V_ind)));
            cr4Ix = intersect(cr1Ix, cr1Ix(ismember(cr1Ix-1,AV_ind)));
        elseif sum(cr1Ix==1)
            fa3Ix = intersect(fa1Ix, fa1Ix(ismember(fa1Ix-1,V_ind)));
            fa4Ix = intersect(fa1Ix, fa1Ix(ismember(fa1Ix-1,AV_ind)));
            cr3Ix = intersect(cr1Ix, cr1Ix(ismember(cr1Ix(2:end)-1,V_ind)));
            cr4Ix = intersect(cr1Ix, cr1Ix(ismember(cr1Ix(2:end)-1,AV_ind)));
        else
            fa3Ix = intersect(fa1Ix, fa1Ix(ismember(fa1Ix-1,V_ind)));
            fa4Ix = intersect(fa1Ix, fa1Ix(ismember(fa1Ix-1,AV_ind)));
            cr3Ix = intersect(cr1Ix, cr1Ix(ismember(cr1Ix-1,V_ind)));
            cr4Ix = intersect(cr1Ix, cr1Ix(ismember(cr1Ix-1,AV_ind)));
        end
        sIxFaIxMatch = [];
        faIxSIxMatch = [];
        mIxCrIxMatch = [];
        faIxMIxMatch = [];
        sIxCrIxMatch = [];
        crIxSIxMatch = [];
        mIxFaIxMatch = [];
        crIxMIxmatch = [];
        faIxCrIxMatch = [];
        crIxFaIxMatch = [];       
        for idir = 1:length(cDirs)
            dirIx = find(tGratingDirectionDeg == cDirs(idir));
            cDirIx = find(tCatchGratingDirectionDeg == cDirs(idir));
            sDir = intersect(Sb1Ix,dirIx);
            mDir = intersect(Mb1Ix,dirIx);
            faDir = intersect(fa1Ix,cDirIx);
            crDir = intersect(cr1Ix,cDirIx);
            nS = length(sDir);
            nM = length(mDir);
            nFA = length(faDir);
            nCR = length(crDir);
            nVec = [nS nM nFA nCR];
            if nS < nFA
                sIxFaIxMatch = [sIxFaIxMatch sDir];
                if any(nVec([1 3]) == 0)
                else
                faIxSIxMatch = [faIxSIxMatch faDir(randperm(nFA,nS))];
                end
            elseif nS >= nFA
                if any(nVec([1 3]) == 0)
                else
                sIxFaIxMatch = [sIxFaIxMatch sDir(randperm(nS,nFA))];
                end
                faIxSIxMatch = [faIxSIxMatch faDir];
            end
            if nM < nFA
                mIxFaIxMatch = [mIxFaIxMatch mDir];
                if any(nVec([2 3]) == 0)
                else
                faIxMIxMatch = [faIxMIxMatch faDir(randperm(nFA,nM))];
                end
            elseif nM >= nFA
                if any(nVec([2 3]) == 0)
                else
                mIxFaIxMatch = [mIxFaIxMatch mDir(randperm(nM,nFA))];
                end
                faIxMIxMatch = [faIxMIxMatch faDir];
            end
            if nS < nCR
                sIxCrIxMatch = [sIxCrIxMatch sDir];
                if any(nVec([1 4]) == 0)
                else
                crIxSIxMatch = [crIxSIxMatch crDir(randperm(nCR,nS))];
                end
            elseif nS >= nCR
                if any(nVec([1 4]) == 0)
                else
                sIxCrIxMatch = [sIxCrIxMatch sDir(randperm(nS,nCR))];
                end
                crIxSIxMatch = [crIxSIxMatch crDir];
            end
            if nM < nCR
                mIxCrIxMatch = [mIxCrIxMatch mDir];
                if any(nVec([2 4]) == 0)
                else
                crIxMIxmatch = [crIxMIxmatch crDir(randperm(nCR,nM))];
                end
            elseif nM >= nCR
                if any(nVec([2 4]) == 0)
                else
                mIxCrIxMatch = [mIxCrIxMatch mDir(randperm(nM,nCR))];
                end
                crIxMIxmatch = [crIxMIxmatch crDir];
            end
            if nCR < nFA
                if any(nVec([3 4]) == 0) 
                else
                faIxCrIxMatch = [faIxCrIxMatch faDir(randperm(nFA,nCR))];
                end
                crIxFaIxMatch = [crIxFaIxMatch crDir];                    
            elseif nCR >= nFA
                if any(nVec([3 4]) == 0)
                else
                crIxFaIxMatch = [crIxFaIxMatch crDir(randperm(nCR,nFA))];
                end
                faIxCrIxMatch = [faIxCrIxMatch faDir];
            end
        end
    end
                
        
        %names of fields
        mouse(imouse).expt(s(:,imouse)).pre_event_frames = pre_event_frames;
        mouse(imouse).expt(s(:,imouse)).post_event_frames = post_event_frames;
        mouse(imouse).expt(s(:,imouse)).visTargets = Dirs;
        if expt(iexp).catch
            mouse(imouse).expt(s(:,imouse)).catchTargets = cDirs;
        end
        mouse(imouse).expt(s(:,imouse)).audTargets = Amps;
        for iDir = 1:length(Dirs)
            mouse(imouse).expt(s(:,imouse)).target(iDir).name = Dirs(iDir);
            mouse(imouse).expt(s(:,imouse)).target(iDir).ind = find(tGratingDirectionDeg==Dirs(iDir));
            if expt(iexp).catch
            mouse(imouse).expt(s(:,imouse)).target(iDir).cind = find(tCatchGratingDirectionDeg==Dirs(iDir));
            end
        end
        
        if isfield(expt,'redChannelOn')
            if expt(iexp).redChannelOn
                mouse(imouse).expt(s(:,imouse)).cells.isLabeledCell = isLabeledCell;
                mouse(imouse).expt(s(:,imouse)).cells.cellType = expt(iexp).redChannelLabel;
            else
                mouse(imouse).expt(s(:,imouse)).cells.isLabeledCell = [];
                mouse(imouse).expt(s(:,imouse)).cells.cellType = [];
            end
        end
        
%         mouse(imouse).expt(s(:,imouse)).cells(3).name = '45';
%         mouse(imouse).expt(s(:,imouse)).cells(4).name = '90';
%         mouse(imouse).expt(s(:,imouse)).cells(5).name = '135';
%         mouse(imouse).expt(s(:,imouse)).cells(6).name = 'untuned';
%         mouse(imouse).expt(s(:,imouse)).cells(7).name = 'all tuned';
%         mouse(imouse).expt(s(:,imouse)).cells(8).name = 'base_excit';
%         mouse(imouse).expt(s(:,imouse)).cells(9).name = 'base_inhib';
%         mouse(imouse).expt(s(:,imouse)).cells(10).name = 'tar_excit';
%         mouse(imouse).expt(s(:,imouse)).cells(11).name = 'tar_inhib';
%         mouse(imouse).expt(s(:,imouse)).cells(12).name = 'base1_resp';
%         mouse(imouse).expt(s(:,imouse)).cells(13).name = 'tar1_resp';
%         mouse(imouse).expt(s(:,imouse)).cells(14).name = 'all cells';
%         mouse(imouse).expt(s(:,imouse)).cells(1).ind = [];
%         mouse(imouse).expt(s(:,imouse)).cells(2).ind = cellsSelect{1};
%         mouse(imouse).expt(s(:,imouse)).cells(3).ind = cellsSelect{2};
%         mouse(imouse).expt(s(:,imouse)).cells(4).ind = cellsSelect{3};
%         mouse(imouse).expt(s(:,imouse)).cells(5).ind = cellsSelect{4};
%         mouse(imouse).expt(s(:,imouse)).cells(6).ind = cellsSelect{5};
%         mouse(imouse).expt(s(:,imouse)).cells(7).ind = cellsSelect{6};
        mouse(imouse).expt(s(:,imouse)).align(1).name = 'press';
        mouse(imouse).expt(s(:,imouse)).align(2).name = 'targetStim';
        mouse(imouse).expt(s(:,imouse)).align(3).name = 'catchStim';
        mouse(imouse).expt(s(:,imouse)).win(1).name = 'pre';
        mouse(imouse).expt(s(:,imouse)).win(2).name = 'trans';
        mouse(imouse).expt(s(:,imouse)).win(1).frames = pre_win;
        mouse(imouse).expt(s(:,imouse)).win(2).frames = trans_win;
        mouse(imouse).expt(s(:,imouse)).info.cyc_time = cycTime;
        mouse(imouse).expt(s(:,imouse)).info.cyc_time_ms = cycTimeMs;
        mouse(imouse).expt(s(:,imouse)).info.nCells = size(dataTC,2);
        mouse(imouse).expt(s(:,imouse)).info.isCatch = expt(iexp).catch;
        mouse(imouse).expt(s(:,imouse)).info.minTrialLengthMs = minTrialLengthMs;
        mouse(imouse).expt(s(:,imouse)).info.minTrialLengthFrames = minTrialLengthFrames;
        mouse(imouse).expt(s(:,imouse)).info.dirs = Dirs;
        mouse(imouse).info.allDirs = Dirs_all;
%         mouse(imouse).expt(s(:,imouse)).tuning(1).name = 'ori preference';
%         mouse(imouse).expt(s(:,imouse)).tuning(2).name = 'OSI';
%         mouse(imouse).expt(s(:,imouse)).tuning(3).name = 'DSI';
%         mouse(imouse).expt(s(:,imouse)).tuning(1).outcome = ori_ind_all;
%         mouse(imouse).expt(s(:,imouse)).tuning(2).outcome = OSI;
%         mouse(imouse).expt(s(:,imouse)).tuning(3).outcome = DSI;
        
        for i = 1:5
            mouse(imouse).expt(s(:,imouse)).align(i).av(1).name = 'visual';
            mouse(imouse).expt(s(:,imouse)).align(i).av(2).name = 'auditory';
            mouse(imouse).expt(s(:,imouse)).align(i).av(3).name = 'visual - prev vis';
            mouse(imouse).expt(s(:,imouse)).align(i).av(4).name = 'visual - prev aud';
            mouse(imouse).expt(s(:,imouse)).align(i).av(5).name = 'auditory - prev aud';
            mouse(imouse).expt(s(:,imouse)).align(i).av(6).name = 'auditory - prev vis';
            for ii = 1:6
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(1).name = 'hit';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(2).name = 'miss';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(3).name = 'FA';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(4).name = 'CR';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(5).name = 'hit_match';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(6).name = 'miss_match';
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(1).n = eval(['length(Sb' num2str(ii) 'Ix)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(2).n = eval(['length(Mb' num2str(ii) 'Ix)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(3).n = eval(['length(Fb' num2str(ii) 'Ix)']);
                if i == 1
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(1).n = eval(['length(SbAR' num2str(ii) 'Ix)']);
                end
                if ii < 3
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(4).n = eval(['length(Rb' num2str(ii) 'Ix)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(5).n = eval(['length(Sb' num2str(ii) 'IxMatch)']);
                mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(6).n = eval(['length(Mb' num2str(ii) 'IxMatch)']);
                end
                if i == 3
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(1).name = 'unmatched hits - dirs';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(2).name = 'unmatched miss - dirs';        
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(3).name = 'miss match FA';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(4).name = 'miss match CR';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(5).name = 'FA match hits';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(6).name = 'FA match miss';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(7).name = 'FA match CR';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(8).name = 'CR match hits';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(9).name = 'CR match miss';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(10).name = 'CR match FA'; 
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(11).name = 'unmatched FA'; 
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(12).name = 'unmatched CR'; 
                end
            end
        end
        
%         %% get correlation and coupling of each neuron to population rate
%         f = mean(dataTC,1);
%         dff = bsxfun(@rdivide, bsxfun(@minus, dataTC, f), f);
%         
%         % Pearson's correleation to population, each neuron
%         timebinS = 0.24;        
%         mouse(imouse).expt(s(:,imouse)).pCorr = popCorr(dff,expt(iexp).frame_rate,timebinS);   
%         
%         % Population coupling
%         mouse(imouse).expt(s(:,imouse)).pCpl = popCouple(dff,expt(iexp).frame_rate);
%         
%         %% get each neuron's fit tuning
%         pref_ind = fitPref(cellFits);
%         mouse(imouse).expt(s(:,imouse)).fitPeak = pref_ind;
%       

        %% Align data to lever down
        ialign = 1;
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataF = zeros(1,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:max(find(cLeverDown+post_event_frames-1 <  size(dataTC,1)),[],2)
            Data(:,:,itrial) = dataTC(cLeverDown(itrial)-pre_event_frames:cLeverDown(itrial)+post_event_frames-1,:);
            DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
        end
        
        %identify trials with motion (large peaks in the derivative)
        if cellsOrDendrites == 1
            motionThreshold = 0.1;
        elseif cellsOrDendrites == 2
            motionThreshold = 0.15;
        elseif cellsOrDendrites == 3
            motionThreshold = 0.15;
        end
        sz = size(DataDFoverF);
        ind_motion = find(max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold);
        mouse(imouse).expt(s(:,imouse)).align(ialign).ind_motion = ind_motion;
        
        figure(motionHist);
        subplot(n,n2,iexp)
        hist(max(diff(squeeze(mean(DataDFoverF,2)),1),[],1));
        hold on
        vline(motionThreshold,'k:')
        xlim([0 0.1])
        xlabel('max diff')
        ylabel('n trials')
               
        %divide data by trial type and outcome
        for iav = 1:2
            if length(eval(['SbAR' num2str(iav) 'Ix']))>0
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).F = Data(:,:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).F = Data(:,:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));        
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).F = Data(:,:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).F = Data(:,:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).F = Data(:,:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)); 
                
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).resp = DataDFoverF(:,:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));        
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = DataDFoverF(:,:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).resp = DataDFoverF(:,:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)); 
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).oneStim = DataDFoverF(:,:,setdiff(eval(['F' num2str(iav) 'oneStim']),ind_motion));                 
                 
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).tcyc = tCyclesOn(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).tcyc = tCyclesOn(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).tcyc = tCyclesOn(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).tcyc = tCyclesOn(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).tcyc = tCyclesOn(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).tcyc = tCyclesOn(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
                
                if iav == 1
                    mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).targetOri = tGratingDirectionDeg(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                    mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).targetOri = tGratingDirectionDeg(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                    mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).targetOri = tGratingDirectionDeg(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                end
                
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).prevTrType = prevTrType(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).prevTrType = prevTrType(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).prevTrType = prevTrType(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).prevTrType = prevTrType(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).prevTrType = prevTrType(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).prevTrType = prevTrType(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));                
                
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).prev2TrType = prev2TrType(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).prev2TrType = prev2TrType(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).prev2TrType = prev2TrType(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).prev2TrType = prev2TrType(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).prev2TrType = prev2TrType(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).prev2TrType = prev2TrType(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
                
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).prevTrOut = prevTrOut(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).prevTrOut = prevTrOut(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).prevTrOut = prevTrOut(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).prevTrOut = prevTrOut(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).prevTrOut = prevTrOut(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).prevTrOut = prevTrOut(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));                
                
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).prev2TrOut = prev2TrOut(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).prev2TrOut = prev2TrOut(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).prev2TrOut = prev2TrOut(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
%                 mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).prev2TrOut = prev2TrOut(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).prev2TrOut = prev2TrOut(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).prev2TrOut = prev2TrOut(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
            else
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).resp = NaN(sz(1), sz(2));
                
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).F = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).F = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).F = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).F = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).F = NaN(sz(1), sz(2));
               
            end
        end
        
        %find trials of a minimum trail length and divide up by outcome
        tCycInd = find(tCyclesOn >= minCyclesAnt);
        for iav = 1:6
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion)));
            if iav < 3
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).cmlvResp = DataDFoverF(1:pre_event_frames+minTrialLengthFrames,:,intersect(tCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
            end
        end
        
        %% find cumulative average for each cycle length
        cycs = 1:max(tCyclesOn);
        for icyc = 1:length(cycs)
            minCyc = cycs(icyc);
            tCycInd = find(tCyclesOn >= minCyc);
            oneCycInd = find(tCyclesOn == minCyc);
            trLengthFrames = cycTime*(icyc);
%         if icyc < 9
        for iav = 1:6
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).cmlvCycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(tCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).cmlvCycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(tCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).cmlvCycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(tCycInd,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion)));
            if iav < 3
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).cmlvCycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(tCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).cmlvCycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(tCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
            end 
        end
        for iav = 1:6
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).cycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(oneCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).cycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).cycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(oneCycInd,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion)));
            if iav < 3
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).cycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(oneCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).cycResp{icyc} = DataDFoverF(1:pre_event_frames+trLengthFrames,:,intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
            end
        end
%         end
        end
        %test for significantly responsive cells
        %to first baseline stimulus (ttest of means)
        [h, p] = ttest(squeeze(mean(DataDFoverF(pre_win,:,setdiff(SIx,ind_motion)),1)), squeeze(mean(DataDFoverF(trans_win,:,setdiff(SIx,ind_motion)),1)), 'dim', 2, 'tail', 'left', 'alpha', 0.05/(length(Dirs)));
        mouse(imouse).expt(s(:,imouse)).align(1).ttest_trans = h;
        %to baseline stimulus before earliest target (ttest of windows across trials)
        [h, p] = ttest(squeeze(mean(DataDFoverF(1:pre_event_frames,:,setdiff(SIx,ind_motion)),1)), squeeze(mean(DataDFoverF(minTime+1:pre_event_frames+minTime,:,setdiff(SIx,ind_motion)),1)), 'dim', 2, 'alpha', 0.05/2);
        mouse(imouse).expt(s(:,imouse)).align(1).ttest_sust = h;
        baseStimRespDiff = squeeze(mean(mean(DataDFoverF(minTime+1:pre_event_frames+minTime,:,setdiff(SIx,ind_motion)),3),1)) - squeeze(mean(mean(DataDFoverF(1:pre_event_frames,:,setdiff(SIx,ind_motion)),3),1));
        %to baseline stimulus before earliest target (vartest of mean response)
        [h, p] = vartest2(squeeze(mean(DataDFoverF(1:pre_event_frames,:,setdiff(SIx,ind_motion)),3)), squeeze(mean(DataDFoverF(minTime+1:pre_event_frames+minTime,:,setdiff(SIx,ind_motion)),3)), 'dim', 2, 'tail', 'left', 'alpha', 0.05/2);
        mouse(imouse).expt(s(:,imouse)).align(1).vartest = h;
        %% Align data to previous stim (could be target or stim before lever release on FA)
        ialign = 2;
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        Data_CR = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDF_CR = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF_CR = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:max(find(cTargetOn+post_event_frames-1 <  size(dataTC,1)),[],2)
            if strcmp(trOut(itrial), 'failure')
                Data(:,:,itrial) = dataTC(cStimOn(itrial)-pre_event_frames:cStimOn(itrial)+post_event_frames-1,:);
                Data_CR(:,:,itrial) = dataTC(cStimOn(itrial)-pre_event_frames-cycTime:cStimOn(itrial)+post_event_frames-cycTime-1,:);
            else
                Data(:,:,itrial) = dataTC(cTargetOn(itrial)-pre_event_frames:cTargetOn(itrial)+post_event_frames-1,:);
                Data_CR(:,:,itrial) = dataTC(cTargetOn(itrial)-pre_event_frames-cycTime:cTargetOn(itrial)+post_event_frames-cycTime-1,:);
            end
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), DataF(:,:,itrial));
            DataDF_CR(:,:,itrial) = bsxfun(@minus, Data_CR(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF_CR(:,:,itrial) = bsxfun(@rdivide, DataDF_CR(:,:,itrial), DataF(:,:,itrial));
        end
        
        %identify trials with motion (large peaks in the derivative)
        ind_motion = find(max(diff(squeeze(mean(DataDFoverF,2)),1),[],1)>motionThreshold);
        mouse(imouse).expt(s(:,imouse)).align(ialign).ind_motion = ind_motion;
        
        %divide data by trial type and outcome
        for iav = 1:6
            if length(eval(['Sb' num2str(iav) 'Ix']))>0
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).resp = DataDFoverF(:,:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = DataDFoverF(:,:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).reactTimes = reactTimes(setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).reactTimes = reactTimes(setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).reactTimes = reactTimes(setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).tcyc = tCyclesOn(:,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).tcyc = tCyclesOn(:,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).tcyc = tCyclesOn(:,setdiff(eval(['Fb' num2str(iav) 'Ix']),ind_motion));
                
                
                if iav < 3
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).resp = DataDFoverF_CR(:,:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).resp = DataDFoverF(:,:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).resp = DataDFoverF(:,:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
                %react times for each trial
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).reactTimes = reactTimes(setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).reactTimes = reactTimes(setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).reactTimes = reactTimes(setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
                %trial length 
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).tcyc = tCyclesOn(:,setdiff(eval(['Rb' num2str(iav) 'Ix']),ind_motion))-1;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).tcyc = tCyclesOn(:,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).tcyc = tCyclesOn(:,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion));
                end
                if expt(iexp).catch && iav == 1
                mouse(imouse).expt(s(:,imouse)).align(3).av(iav).outcome(1).resp = DataDFoverF(:,:,setdiff(sIxFaIxMatch,ind_motion)); %'hits match FA';
                mouse(imouse).expt(s(:,imouse)).align(3).av(iav).outcome(2).resp = DataDFoverF(:,:,setdiff(sIxCrIxMatch,ind_motion)); %'hits match CR';        
                mouse(imouse).expt(s(:,imouse)).align(3).av(iav).outcome(3).resp = DataDFoverF(:,:,setdiff(mIxFaIxMatch,ind_motion)); %'miss match FA';
                mouse(imouse).expt(s(:,imouse)).align(3).av(iav).outcome(4).resp = DataDFoverF(:,:,setdiff(mIxCrIxMatch,ind_motion)); %'miss match CR';
                end
            else
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).resp = NaN(sz(1), sz(2));
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).tcyc = NaN;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).tcyc = NaN;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).tcyc = NaN;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).tcyc = NaN;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).tcyc = NaN;
                mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).tcyc = NaN;
                end
        end
        
        %divide data by target direction
        mouse(imouse).expt(s(:,imouse)).align(ialign).ttest_trans = zeros(size(dataTC,2), length(Dirs));
        tarStimRespDiff = zeros(size(dataTC,2), length(Dirs));
        for iDir = 1:length(Dirs)
            for iav = [1 3 4];
            vInd = setdiff(find(tGratingDirectionDeg==Dirs(iDir)), ind_motion);
            %if iDir == 1, then the direction is 0 and it's an auditory trial
            if iDir == 1
                indS = intersect(vInd, Sb2Ix);
                indM = intersect(vInd, Mb2Ix);
                if iav < 3
                indSM = intersect(vInd, Sb2IxMatch);
                indMM = intersect(vInd, Mb2IxMatch);
                end
            else
                indS = intersect(vInd, eval(['Sb' num2str(iav) 'Ix']));%Sb1Ix);
                indM = intersect(vInd, eval(['Mb' num2str(iav) 'Ix']));%Mb1Ix);
                if iav < 3
                indSM = intersect(vInd, Sb1IxMatch);
                indMM = intersect(vInd, Mb1IxMatch);
                end
            end
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).stimResp{iDir} = DataDFoverF(:,:,indS);
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).stimResp{iDir} = DataDFoverF(:,:,indM);  
            if iav < 3
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).stimResp{iDir} = DataDFoverF(:,:,indSM);
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).stimResp{iDir} = DataDFoverF(:,:,indMM);
            end
            
             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).trialL{iDir} = tCyclesOn(:,indS);
             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).trialL{iDir} = tCyclesOn(:,indM);
             
            %test for responsiveness to any one of the presented directions
            if iav == 1
            [h, p] = ttest(squeeze(mean(DataDFoverF(pre_win,:,[indS indM]),1)), squeeze(mean(DataDFoverF(trans_win,:,[indS indM]),1)), 'dim', 2, 'tail', 'left', 'alpha', 0.05/(length(Dirs)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).ttest_trans(:,iDir) = h;
            tarStimRespDiff(:,iDir) = squeeze(mean(mean(DataDFoverF(trans_win,:,[indS indM]),3),1)) - squeeze(mean(mean(DataDFoverF(pre_win,:,[indS indM]),3),1));
            end
            end
        end
        
        %divide data by target sound amplitude
        for iAmp = 1:length(Amps)
            iav = 2;
            if isfield(input, 'tSoundTargetAmplitude')
                aInd = setdiff(find(tSoundTargetAmplitude==Amps(iAmp)), ind_motion)';
            else
                aInd = setdiff(find(celleqel2mat_padded(input.tDoAuditoryDetect(tr))), ind_motion);
            end
            %if iDir == 1, then the amplitude is 0 and it's a visual trial
            indS = NaN; indM = NaN; indSM = NaN; indMM = NaN;
            if iAmp == 1
                indS = intersect(aInd, Sb1Ix);
                indM = intersect(aInd, Mb1Ix);
                indSM = intersect(aInd, Sb1IxMatch);
                indMM = intersect(aInd, Mb1IxMatch);
            else
                indS = intersect(aInd, Sb2Ix);
                indM = intersect(aInd, Mb2Ix);
                indSM = intersect(aInd, Sb2IxMatch);
                indMM = intersect(aInd, Mb2IxMatch);
            end
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).stimResp{iAmp} = DataDFoverF(:,:,indS);
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).stimResp{iAmp} = DataDFoverF(:,:,indM);        
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).stimResp{iAmp} = DataDFoverF(:,:,indSM);
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).stimResp{iAmp} = DataDFoverF(:,:,indMM);
        end
        
        %% Align data to vis catch stim
        if expt(iexp).catch
        ialign = 3;
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:max(find(cCatchOn+post_event_frames-1 <  size(dataTC,1)),[],2)
            if expt(iexp).catch
            if strcmp(catchTrialOutcome{itrial}, 'failure') | isempty(catchTrialOutcome{itrial})
                Data(:,:,itrial) = NaN(pre_event_frames+post_event_frames,size(dataTC,2));
            else
                Data(:,:,itrial) = dataTC(cCatchOn(itrial)-pre_event_frames:cCatchOn(itrial)+post_event_frames-1,:);
            end
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), DataF(:,:,itrial));
            end
        end
        
%         %identify trials with motion (large peaks in the derivative)
%         ind_motion = find(max(diff(squeeze(nanmean(DataDFoverF,2)),1),[],1)>motionThreshold);
%         disp(length(ind_motion));
%         mouse(imouse).expt(s(:,imouse)).align(ialign).ind_motion = ind_motion;
                
        %divide data by trial type and outcome
        iav = 1;
   
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).resp = DataDFoverF(:,:,setdiff(faIxSIxMatch,ind_motion)); %'FA match hits';
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).resp = DataDFoverF(:,:,setdiff(faIxMIxMatch,ind_motion)); %'FA match miss';
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(7).resp = DataDFoverF(:,:,setdiff(faIxCrIxMatch,ind_motion)); %'FA match CR';
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(8).resp = DataDFoverF(:,:,setdiff(crIxSIxMatch,ind_motion)); %'CR match hits';
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(9).resp = DataDFoverF(:,:,setdiff(crIxMIxmatch,ind_motion)); %'CR match miss';
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(10).resp = DataDFoverF(:,:,setdiff(crIxFaIxMatch,ind_motion)); %'CR match FA';
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(11).resp = DataDFoverF(:,:,setdiff(fa1Ix,ind_motion)); %'unmatched FA';        
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(12).resp = DataDFoverF(:,:,setdiff(cr1Ix,ind_motion)); %'unmatched CR';
        
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(11).tcyc = catchCyclesOn(:,setdiff(fa1Ix,ind_motion)); %'unmatched FA';        
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(12).tcyc = catchCyclesOn(:,setdiff(cr1Ix,ind_motion)); %'unmatched CR';
        %divide data by catch direction
        mouse(imouse).expt(s(:,imouse)).info.cDirs = cDirs;
        for iav = [1 3 4]
        cind = [];
        for iDir = 1:length(cDirs)
            cind = find(tCatchGratingDirectionDeg == cDirs(iDir));
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).stimResp{iDir} = mouse(imouse).expt(s(:,imouse)).align(2).av(iav).outcome(1).stimResp{iDir};
%             mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).stimResp{iDir} = mouse(imouse).expt(s(:,imouse)).align(2).av(iav).outcome(2).stimResp{iDir};
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).stimResp{iDir} = DataDFoverF(:,:,setdiff(intersect(eval(['fa' num2str(iav) 'Ix']),cind),ind_motion));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).stimResp{iDir} = DataDFoverF(:,:,setdiff(intersect(eval(['cr' num2str(iav) 'Ix']),cind),ind_motion));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).tcyc{iDir} = catchCyclesOn(:,setdiff(intersect(eval(['fa' num2str(iav) 'Ix']),cind),ind_motion));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).tcyc{iDir} = catchCyclesOn(:,setdiff(intersect(eval(['cr' num2str(iav) 'Ix']),cind),ind_motion));
            
        end
        
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).FAreacttime = catchTrialReactTimeMs(eval(['fa' num2str(iav) 'Ix']));
        
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).stimResp = mouse(imouse).expt(s(:,imouse)).align(2).av(iav).outcome(1).stimResp(find(ismember(Dirs,cDirs)));
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).stimResp = mouse(imouse).expt(s(:,imouse)).align(2).av(iav).outcome(2).stimResp(find(ismember(Dirs,cDirs)));
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).tcyc = mouse(imouse).expt(s(:,imouse)).align(2).av(iav).outcome(1).trialL(find(ismember(Dirs,cDirs)));
        mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).tcyc = mouse(imouse).expt(s(:,imouse)).align(2).av(iav).outcome(2).trialL(find(ismember(Dirs,cDirs)));
        
        end
%         for idir = 1:length(cDirs)
%            dirIx = find(tCatchGratingDirectionDeg == cDirs(idir));
%            faIxDir = intersect(fa1Ix,dirIx);
%            crIxDir = intersect(cr1Ix,dirIx);
%            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).stimResp{idir} = DataDFoverF(:,:,faIxDir);
%            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).stimResp{idir} = DataDFoverF(:,:,crIxDir);
%         end        
        else
            mouse(imouse).expt(s(:,imouse)).info.cDirs = [];
        end
        
        %% Align data to lever down, include target
        ialign = 4;
        Data = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataF = zeros(1,size(dataTC,2),ntrials);
        DataDF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        DataDFoverF = zeros(pre_event_frames+post_event_frames,size(dataTC,2),ntrials);
        for itrial = 1:max(find(cLeverDown+post_event_frames-1 <  size(dataTC,1)),[],2)
            Data(:,:,itrial) = dataTC(cLeverDown(itrial)-pre_event_frames:cLeverDown(itrial)+post_event_frames-1,:);
            DataF(:,:,itrial) = mean(Data(1:pre_event_frames,:,itrial),1);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), DataF(:,:,itrial));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:pre_event_frames,:,itrial),1));
        end
        
        cycs = 1:max(tCyclesOn);
        for icyc = 1:length(cycs)
            minCyc = cycs(icyc);
            oneCycInd = find(tCyclesOn == minCyc);
        for iav = 1:2
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).cycResp{icyc} = DataDFoverF(:,:,intersect(oneCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).cycResp{icyc} = DataDFoverF(:,:,intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));              
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).cycResp{icyc} = DataDFoverF(:,:,intersect(oneCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).cycResp{icyc} = DataDFoverF(:,:,intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
            
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).tTargetOn{icyc} = cTargetOn(intersect(oneCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).tTargetOn{icyc} = cTargetOn(intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).tTargetOn{icyc} = cTargetOn(intersect(oneCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).tTargetOn{icyc} = cTargetOn(intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
            
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).direction{icyc} = tGratingDirectionDeg(intersect(oneCycInd,setdiff(eval(['SbAR' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).direction{icyc} = tGratingDirectionDeg(intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).direction{icyc} = tGratingDirectionDeg(intersect(oneCycInd,setdiff(eval(['Sb' num2str(iav) 'IxMatch']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).direction{icyc} = tGratingDirectionDeg(intersect(oneCycInd,setdiff(eval(['Mb' num2str(iav) 'IxMatch']),ind_motion)));
                        
            if expt(iexp).catch & iav == 1
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).cycResp{icyc} = DataDFoverF(:,:,intersect(oneCycInd,setdiff(eval(['fa' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).cycResp{icyc} = DataDFoverF(:,:,intersect(oneCycInd,setdiff(eval(['cr' num2str(iav) 'Ix']),ind_motion))); 
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).tTargetOn{icyc} = cCatchOn(intersect(oneCycInd,setdiff(eval(['fa' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).tTargetOn{icyc} = cCatchOn(intersect(oneCycInd,setdiff(eval(['cr' num2str(iav) 'Ix']),ind_motion)));  
            
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).direction{icyc} = tCatchGratingDirectionDeg(intersect(oneCycInd,setdiff(eval(['fa' num2str(iav) 'Ix']),ind_motion)));
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).direction{icyc} = tCatchGratingDirectionDeg(intersect(oneCycInd,setdiff(eval(['cr' num2str(iav) 'Ix']),ind_motion)));  
            
            end
            
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(1).name = 'hit';
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(2).name = 'miss'; 
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(3).name = 'invalid hit';
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(4).name = 'invalid miss';              
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(5).name = 'hit - matched';
            mouse(imouse).expt(s(:,imouse)).align(ialign).av(iav).outcome(6).name = 'miss - matched';
            
        end
        end
%         %% find cells that are responsive to at least one of the target
%         %directions or the first baseline visual stimulus 
%         base_ind = find(mouse(imouse).expt(s(:,imouse)).align(1).ttest_trans);
%         resp_ind = find(sum(mouse(imouse).expt(s(:,imouse)).align(2).ttest_trans,2)>=1);
%         mouse(imouse).expt(s(:,imouse)).cells(1).ind = unique([base_ind; resp_ind]);
%         mouse(imouse).expt(s(:,imouse)).cells(8).ind = intersect(find(baseStimRespDiff > 0),(find(mouse(imouse).expt(s(:,imouse)).align(1).ttest_sust)));
%         mouse(imouse).expt(s(:,imouse)).cells(9).ind = intersect(find(baseStimRespDiff < 0),(find(mouse(imouse).expt(s(:,imouse)).align(1).ttest_sust)));
%         mouse(imouse).expt(s(:,imouse)).cells(10).ind = intersect(find(mean(tarStimRespDiff,2) > 0),resp_ind);
%         mouse(imouse).expt(s(:,imouse)).cells(11).ind = intersect(find(mean(tarStimRespDiff,2) < 0),resp_ind);
%         mouse(imouse).expt(s(:,imouse)).cells(12).ind = base_ind;
%         mouse(imouse).expt(s(:,imouse)).cells(13).ind = resp_ind;
%         if isempty(mouse(imouse).expt(s(:,imouse)).cells(1).ind) & isempty(mouse(imouse).expt(s(:,imouse)).cells(8).ind)
%             mouse(imouse).expt(s(:,imouse)).cells(14).ind = mouse(imouse).expt(s(:,imouse)).cells(1).ind;
%         elseif isempty(mouse(imouse).expt(s(:,imouse)).cells(1).ind)
%             mouse(imouse).expt(s(:,imouse)).cells(14).ind = mouse(imouse).expt(s(:,imouse)).cells(8).ind;
%         elseif isempty(mouse(imouse).expt(s(:,imouse)).cells(8).ind)
%             mouse(imouse).expt(s(:,imouse)).cells(14).ind = mouse(imouse).expt(s(:,imouse)).cells(1).ind;
%         else
%             mouse(imouse).expt(s(:,imouse)).cells(14).ind = unique(cat(1,mouse(imouse).expt(s(:,imouse)).cells(1).ind,mouse(imouse).expt(s(:,imouse)).cells(8).ind));
%         end
    end
    
    figure(motionHist)
    if cellsOrDendrites == 1
        if doRawData == 1
            save(fullfile(rc.caOutputDir, dataGroup, [mouse_str '_rawDataTC_CaSummary_cells' datasetStr(5:end) '.mat']), 'mouse','-v7.3');
        else
            save(fullfile(rc.caOutputDir, dataGroup, [mouse_str '_CaSummary_cells' datasetStr(5:end) '.mat']), 'mouse','-v7.3');
            print(fullfile(rc.caOutputDir, dataGroup, [datasetStr(5:end) '_motionHist.pdf']), '-dpdf')
        end
    elseif cellsOrDendrites == 2
        if doRawData == 1
            save(fullfile(rc.caOutputDir, dataGroup, [mouse_str '_rawDataTC_CaSummary_dendrites' datasetStr(5:end) '.mat']), 'mouse','-v7.3');
        else
            save(fullfile(rc.caOutputDir, dataGroup, [mouse_str '_CaSummary_dendrites' datasetStr(5:end) '.mat']), 'mouse','-v7.3');
            print(fullfile(rc.caOutputDir, dataGroup, [datasetStr(5:end) '_motionHist.pdf']), '-dpdf')
        end
    else
        print(fullfile(rc.caOutputDir, dataGroup, [datasetStr(5:end) '_motionHist.pdf']), '-dpdf')
        save(fullfile(rc.caOutputDir, dataGroup, [mouse_str '_CaSummary' datasetStr(5:end) '.mat']), 'mouse','-v7.3');
    end
end
        
        
        

