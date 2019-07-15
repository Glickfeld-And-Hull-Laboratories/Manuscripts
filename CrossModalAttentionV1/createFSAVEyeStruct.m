clear all
close all
ds = 'FSAV_attentionV1';
eval(ds)
rc = behavConstsAV;
imgParams_FSAV

titleStr = ds(6:end);
fnout = fullfile(rc.caOutputDir, ds,[titleStr '_eye_']);
%%
mice = unique({expt.SubNum});
nMice = length(mice);
nExpt = size(expt,2);

preEventFr = round(preEventMs_eye*frameRateHz/1000);
postEventFr = round(postEventMs_eye*frameRateHz/1000);

%%
msEye = struct;

for iexp = 1:nExpt
    if iexp == 1
        nExptPerMs = zeros(1,nMice);
    end
        ms = expt(iexp).mouse;
        subnum = expt(iexp).SubNum;
        im = find(strcmp(mice,subnum));
        msEye(im).name = mice{im};
        nExptPerMs(im) = nExptPerMs(im)+1;
    if ~isnan(expt(iexp).eyeradrange)
        if strcmp(ds,'FSAV_attentionV1')
            msEye(im).eyeSameAsImgDay = true;
        end
        dt = expt(iexp).date;
        fprintf('Expt#%s: %s-%s\n',num2str(iexp),ms,dt)
        msEye(im).expt(nExptPerMs(im)).date = dt;
        
        runFolders = expt(iexp).runs;
        runTimes = expt(iexp).time_mat;
        nRun = size(runFolders,1);
        fprintf('Concatenating runs...\n')
        for irun = 1:nRun
            if irun == 1
                nFramesPerRun = zeros(1,nRun);
                AreaAll = [];
                CentroidAll = [];      
                mw = [];
            end
            fn = fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,...
                runFolders(irun,:));
            load(fullfile(fn,'eyeTC.mat'))
            if ~isempty(expt(iexp).nFramesPerRun)
                if ~isnan(expt(iexp).nFramesPerRun{irun})
                    nFr = expt(iexp).nFramesPerRun{irun};
                else
                load(fullfile(rc.ashleyData,ms,'two-photon imaging',dt,...
                    runFolders(irun,:),[runFolders(irun,:) '_000_000.mat']))
                nFr = info.config.frames;
                end
            else
                load(fullfile(rc.ashleyData,ms,'two-photon imaging',dt,...
                    runFolders(irun,:),[runFolders(irun,:) '_000_000.mat']))
                nFr = info.config.frames;
            end
            nFramesPerRun(irun) = nFr;
            AreaAll = cat(1,AreaAll,Area(1:nFr,:));
            CentroidAll = cat(1,CentroidAll,Centroid(1:nFr,:));
            
            mw = [mw loadMworksFile(subnum,dt,runTimes(irun,:))];
        end
        mw = concatenateDataBlocks(mw);
        if isnan(expt(iexp).trial_range)
            tr = 1:length(mw.trialOutcomeCell);
        else
            tr = expt(iexp).trial_range;
        end    
        clear Radius Area Centroid
        radius = sqrt(AreaAll./pi).*eyeMmPerPix;
        horzPos = CentroidAll(:,1).*eyeMmPerPix;
        vertPos = CentroidAll(:,2).*eyeMmPerPix;
        nfr = length(radius);
        
        fprintf('Aligning data...\n')
        %behavior info needed
        cumTrialsPerRun = cumsum(mw.trialsSinceReset);
        for irun = 1:nRun
            if irun == 1
                offset = 0;
                tLeverDown = celleqel2mat_padded(mw.cLeverDown);
                tTargetOn = celleqel2mat_padded(mw.cTargetOn);
            else
                offset = offset+nFramesPerRun(irun-1);
                ind = (cumTrialsPerRun(irun-1)+1):cumTrialsPerRun(irun);
                tLeverDown(ind) = tLeverDown(ind)+offset;
                tTargetOn(ind) = tTargetOn(ind)+offset;
            end
        end
        tLeverDown = tLeverDown(tr);
        tTargetOn = tTargetOn(tr);        
        tCyclesOn = cell2mat(mw.tCyclesOn(tr));
        trialOutcome = mw.trialOutcomeCell(tr);
        if expt(iexp).catch
            isCatchTrial = ~isnan(celleqel2mat_padded(mw.cCatchOn(tr)));
        end
        isVisTrial = celleqel2mat_padded(mw.tGratingDirectionDeg(tr)) > 0;
        
        %behavior info just for eye data
        if iscell(mw.nFramesOn)
            on = unique(cell2mat(mw.nFramesOn));
            off = unique(cell2mat(mw.nFramesOff));
            cycTimeFr = on+off;
        else
%             cycTimeFr = mw.nFramesOn + mw.nFramesOff;
            cycTimeFr = round(((mw.stimOnTimeMs+mw.stimOffTimeMs)...
                ./1000)*frameRateHz);
        end
        hits = strcmp(trialOutcome,'success') & ~isCatchTrial;
        nAllTrials = length(hits);
        nHitTrials = sum(hits);
        trialType = nan(1,nHitTrials);
        trialType(isVisTrial(hits)) = visualTrials;
        trialType(~isVisTrial(hits)) = auditoryTrials;
        
        msEye(im).expt(nExptPerMs(im)).pos(1).name = 'radius';
        msEye(im).expt(nExptPerMs(im)).pos(2).name = 'horizontal position';
        msEye(im).expt(nExptPerMs(im)).pos(3).name = 'vertical position';
        
        %start align
        rad_align = nan(preEventFr+postEventFr,nHitTrials);
        horzPos_align = nan(preEventFr+postEventFr,nHitTrials);
        vertPos_align = nan(preEventFr+postEventFr,nHitTrials);
        for i = 1:nAllTrials
            if i == 1
                itrial = 0;
            end
            if hits(i)
                itrial = itrial+1;
                ind = (tLeverDown(i) - preEventFr):(tLeverDown(i)+postEventFr-1);
                rad_align(:,itrial) = radius(ind);
                horzPos_align(:,itrial) = horzPos(ind);
                vertPos_align(:,itrial) = vertPos(ind);
            end
        end
        
        rad_align_norm = rad_align./mean(rad_align(1:preEventFr,:),1);
        horzPos_align_sub = horzPos_align - ...
            mean(horzPos_align(1:preEventFr,:),1);
        vertPos_align_sub = vertPos_align - ...
            mean(vertPos_align(1:preEventFr,:),1);
        
        shortTrialFr = shortTrialTimeS_eye*frameRateHz;
        minCyc = ceil(shortTrialFr/cycTimeFr);
        minCycInd = tCyclesOn(hits) >= minCyc;
        
        for iav = 1:2
            msEye(im).expt(nExptPerMs(im)).pos(1).av(iav)...
                .align(alignStart).shortTC = ...
                rad_align_norm(1:(preEventFr+shortTrialFr),...
                minCycInd & trialType == iav);
            msEye(im).expt(nExptPerMs(im)).pos(2).av(iav)...
                .align(alignStart).shortTC = ...
                horzPos_align_sub(1:(preEventFr+shortTrialFr),...
                minCycInd & trialType == iav);
            msEye(im).expt(nExptPerMs(im)).pos(3).av(iav)...
                .align(alignStart).shortTC = ...
                vertPos_align_sub(1:(preEventFr+shortTrialFr),...
                minCycInd & trialType == iav);
        end
        
        longTrialFr = longTrialTimeS_eye*frameRateHz;
        minCyc = ceil(longTrialFr/cycTimeFr);
        minCycInd = tCyclesOn(hits) >= minCyc;
        
        for iav = 1:2
            msEye(im).expt(nExptPerMs(im)).pos(1).av(iav)...
                .align(alignStart).longTC = ...
                rad_align_norm(1:(preEventFr+longTrialFr),...
                minCycInd & trialType == iav);
            msEye(im).expt(nExptPerMs(im)).pos(2).av(iav)...
                .align(alignStart).longTC = ...
                horzPos_align_sub(1:(preEventFr+longTrialFr),...
                minCycInd & trialType == iav);
            msEye(im).expt(nExptPerMs(im)).pos(3).av(iav)...
                .align(alignStart).longTC = ...
                vertPos_align_sub(1:(preEventFr+longTrialFr),...
                minCycInd & trialType == iav);
        end
        
        %target align
        targetFr = targetTimeS_eye*frameRateHz;
        
        rad_align = nan(preEventFr+targetFr,nHitTrials);
        horzPos_align = nan(preEventFr+targetFr,nHitTrials);
        vertPos_align = nan(preEventFr+targetFr,nHitTrials);
        for i = 1:nAllTrials
            if i == 1
                itrial = 0;
            end
            if hits(i)
                itrial = itrial+1;
                ind = (tTargetOn(i) - preEventFr):(tTargetOn(i)+targetFr-1);
                rad_align(:,itrial) = radius(ind);
                horzPos_align(:,itrial) = horzPos(ind);
                vertPos_align(:,itrial) = vertPos(ind);
            end
        end
        
        rad_align_norm = rad_align./mean(rad_align(1:preEventFr,:),1);
        horzPos_align_sub = horzPos_align - ...
            mean(horzPos_align(1:preEventFr,:),1);
        vertPos_align_sub = vertPos_align - ...
            mean(vertPos_align(1:preEventFr,:),1);
                
        for iav = 1:2
            msEye(im).expt(nExptPerMs(im)).pos(1).av(iav)...
                .align(alignTarget).TC = ...
                rad_align_norm(1:(preEventFr+targetFr),trialType == iav);
            msEye(im).expt(nExptPerMs(im)).pos(2).av(iav)...
                .align(alignTarget).TC = ...
                horzPos_align_sub(1:(preEventFr+targetFr),trialType == iav);
            msEye(im).expt(nExptPerMs(im)).pos(3).av(iav)...
                .align(alignTarget).TC = ...
                vertPos_align_sub(1:(preEventFr+targetFr),trialType == iav);
        end
        
    else
        if strcmp(ds,'FSAV_attentionV1')
            AWEyeDatasets_AW
            msEye(im).eyeSameAsImgDay = false;
            if nExptPerMs(im) == 1
                fprintf('Using eye-only experiments for mouse %s\n',ms)
                eyeExptInd = strcmp({eyeExpt.SubNum},mice(im));
                eyeExptMs = eyeExpt(eyeExptInd);
                nExptMs = sum(eyeExptInd);
                for imsexpt = 1:nExptMs    
                    dt = eyeExptMs(imsexpt).date;
                    fprintf('Eye-Only Expt#%s: %s-%s\n',num2str(imsexpt),ms,dt)
                    msEye(im).expt(imsexpt).date = dt;
                    
                    fprintf('Concatenating runs...\n')
                    runFolders = eyeExptMs(imsexpt).runs;
                    runTimes = eyeExptMs(imsexpt).time_mat;
                    nRun = size(runFolders,1);
                    fn = fullfile(rc.ashleyAnalysis,ms,'eye tracking',dt);                    
                    for irun = 1:nRun
                        if irun == 1
                            eyeFileName = sprintf('%s-%s',ms,dt);
                        end
                        if irun == nRun
                            eyeFileName = sprintf('%s-%s_pupil.mat',eyeFileName,...
                                runFolders(irun,:));
                        else
                            eyeFileName = sprintf('%s-%s',eyeFileName,...
                                runFolders(irun,:));
                        end
                    end
                    load(fullfile(fn,eyeFileName))
                    for irun = 1:nRun
                        if irun == 1
                            AreaAll = [];
                            CentroidAll = [];      
                            mw = [];
                        end
                        load(fullfile(rc.ashleyData,ms,'eye tracking',dt,...
                            runFolders(irun,:),[runFolders(irun,:) '_000_000.mat']))
                        if ~isempty(eyeExptMs(imsexpt).nFramesPerRun)
                            if ~isnan(eyeExptMs(imsexpt).nFramesPerRun{irun})
                                nFr = eyeExptMs(imsexpt).nFramesPerRun{irun};
                            else
                            load(fullfile(rc.ashleyData,ms,'eye tracking',dt,...
                                runFolders(irun,:),[runFolders(irun,:) '_000_000.mat']))
                            nFr = info.config.frames;
                            end
                        else
                            load(fullfile(rc.ashleyData,ms,'eye tracking',dt,...
                                runFolders(irun,:),[runFolders(irun,:) '_000_000.mat']))
                            nFr = info.config.frames;
                        end
                        AreaAll = cat(1,AreaAll,Area{irun}(1:nFr,:));
                        CentroidAll = cat(1,CentroidAll,Centroid{irun}(1:nFr,:));
                        mw = [mw loadMworksFile(subnum,dt,runTimes(irun,:))];
                    end
                    mw = concatenateDataBlocks(mw);
                    
                    if isnan(eyeExptMs(imsexpt).trial_range)
                        tr = 1:length(mw.trialOutcomeCell);
                    else
                        tr = eyeExptMs(imsexpt).trial_range;
                    end    
                    clear Radius Area Centroid
                    radius = sqrt(AreaAll./pi).*eyeMmPerPix;
                    horzPos = CentroidAll(:,1).*eyeMmPerPix;
                    vertPos = CentroidAll(:,2).*eyeMmPerPix;
                    nfr = length(radius);
                    
                    fprintf('Aligning data...\n')
                    %behavior info needed
                    cumTrialsPerRun = cumsum(mw.trialsSinceReset);
                    for irun = 1:nRun
                        if irun == 1
                            offset = 0;
                            tLeverDown = celleqel2mat_padded(mw.cLeverDown);
                            tTargetOn = celleqel2mat_padded(mw.cTargetOn);
                        else
                            offset = offset+cumTrialsPerRun(irun-1);
                            ind = (cumTrialsPerRun(irun-1)+1):cumTrialsPerRun(irun);
                            tLeverDown(ind) = tLeverDown(ind)+offset;
                            tTargetOn(ind) = tTargetOn(ind)+offset;
                        end
                    end
                    tLeverDown = tLeverDown(tr);
                    tTargetOn = tTargetOn(tr);
                    tCyclesOn = cell2mat(mw.tCyclesOn(tr));
                    trialOutcome = mw.trialOutcomeCell(tr);
                    isCatchTrial = ~isnan(celleqel2mat_padded(mw.cCatchOn(tr)));
                    isVisTrial = celleqel2mat_padded(mw.tGratingDirectionDeg(tr)) > 0;
                    

                    %behavior info just for eye data
                    if iscell(mw.nFramesOn)
                        on = unique(cell2mat(mw.nFramesOn));
                        off = unique(cell2mat(mw.nFramesOff));
                        cycTimeFr = on+off;
                    else
%                         cycTimeFr = mw.nFramesOn + mw.nFramesOff;
                        cycTimeFr = round(((mw.stimOnTimeMs+mw.stimOffTimeMs)...
                            ./1000)*frameRateHz);
                    end
                    hits = strcmp(trialOutcome,'success') & ~isCatchTrial;
                    nAllTrials = length(hits);
                    nHitTrials = sum(hits);
                    trialType = nan(1,nHitTrials);
                    trialType(isVisTrial(hits)) = visualTrials;
                    trialType(~isVisTrial(hits)) = auditoryTrials;

                    msEye(im).expt(imsexpt).pos(1).name = 'radius';
                    msEye(im).expt(imsexpt).pos(2).name = 'horizontal position';
                    msEye(im).expt(imsexpt).pos(3).name = 'vertical position';
                    
                    %start align
                    rad_align = nan(preEventFr+postEventFr,nHitTrials);
                    horzPos_align = nan(preEventFr+postEventFr,nHitTrials);
                    vertPos_align = nan(preEventFr+postEventFr,nHitTrials);
                    for i = 1:nAllTrials
                        if i == 1
                            itrial = 0;
                        end
                        if hits(i)
                            itrial = itrial+1;
                            ind = (tLeverDown(i) - preEventFr):(tLeverDown(i)+postEventFr-1);
                            rad_align(:,itrial) = radius(ind);
                            horzPos_align(:,itrial) = horzPos(ind);
                            vertPos_align(:,itrial) = vertPos(ind);
                        end
                    end

                    rad_align_norm = rad_align./mean(rad_align(1:preEventFr,:),1);
                    horzPos_align_sub = horzPos_align - ...
                        mean(horzPos_align(1:preEventFr,:),1);
                    vertPos_align_sub = vertPos_align - ...
                        mean(vertPos_align(1:preEventFr,:),1);

                    shortTrialFr = shortTrialTimeS_eye*frameRateHz;
                    minCyc = ceil(shortTrialFr/cycTimeFr);
                    minCycInd = tCyclesOn(hits) >= minCyc;

                    for iav = 1:2
                        msEye(im).expt(imsexpt).pos(1).av(iav)...
                            .align(alignStart).shortTC = ...
                            rad_align_norm(1:(preEventFr+shortTrialFr),...
                            minCycInd & trialType == iav);
                        msEye(im).expt(imsexpt).pos(2).av(iav)...
                            .align(alignStart).shortTC = ...
                            horzPos_align_sub(1:(preEventFr+shortTrialFr),...
                            minCycInd & trialType == iav);
                        msEye(im).expt(imsexpt).pos(3).av(iav)...
                            .align(alignStart).shortTC = ...
                            vertPos_align_sub(1:(preEventFr+shortTrialFr),...
                            minCycInd & trialType == iav);
                    end

                    longTrialFr = longTrialTimeS_eye*frameRateHz;
                    minCyc = ceil(longTrialFr/cycTimeFr);
                    minCycInd = tCyclesOn(hits) >= minCyc;

                    for iav = 1:2
                        msEye(im).expt(imsexpt).pos(1).av(iav)...
                            .align(alignStart).longTC = ...
                            rad_align_norm(1:(preEventFr+longTrialFr),...
                            minCycInd & trialType == iav);
                        msEye(im).expt(imsexpt).pos(2).av(iav)...
                            .align(alignStart).longTC = ...
                            horzPos_align_sub(1:(preEventFr+longTrialFr),...
                            minCycInd & trialType == iav);
                        msEye(im).expt(imsexpt).pos(3).av(iav)...
                            .align(alignStart).longTC = ...
                            vertPos_align_sub(1:(preEventFr+longTrialFr),...
                            minCycInd & trialType == iav);
                    end
                    
                    %target align
                    targetFr = targetTimeS_eye*frameRateHz;

                    rad_align = nan(preEventFr+targetFr,nHitTrials);
                    horzPos_align = nan(preEventFr+targetFr,nHitTrials);
                    vertPos_align = nan(preEventFr+targetFr,nHitTrials);
                    for i = 1:nAllTrials
                        if i == 1
                            itrial = 0;
                        end
                        if hits(i)
                            itrial = itrial+1;
                            ind = (tTargetOn(i) - preEventFr):(tTargetOn(i)+targetFr-1);
                            rad_align(:,itrial) = radius(ind);
                            horzPos_align(:,itrial) = horzPos(ind);
                            vertPos_align(:,itrial) = vertPos(ind);
                        end
                    end

                    rad_align_norm = rad_align./mean(rad_align(1:preEventFr,:),1);
                    horzPos_align_sub = horzPos_align - ...
                        mean(horzPos_align(1:preEventFr,:),1);
                    vertPos_align_sub = vertPos_align - ...
                        mean(vertPos_align(1:preEventFr,:),1);

                    for iav = 1:2
                        msEye(im).expt(imsexpt).pos(1).av(iav)...
                            .align(alignTarget).TC = ...
                            rad_align_norm(1:(preEventFr+targetFr),trialType == iav);
                        msEye(im).expt(imsexpt).pos(2).av(iav)...
                            .align(alignTarget).TC = ...
                            horzPos_align_sub(1:(preEventFr+targetFr),trialType == iav);
                        msEye(im).expt(imsexpt).pos(3).av(iav)...
                            .align(alignTarget).TC = ...
                            vertPos_align_sub(1:(preEventFr+targetFr),trialType == iav);
                    end
                    msEye(im).expt(imsexpt).pos(1).av(iav)...
                            .align(alignStart).name = 'align to start';
                    msEye(im).expt(imsexpt).pos(1).av(iav)...
                            .align(alignTarget).name = 'align to target';
                end
            end
        else
            fprintf('Add pupil datasets for %s\n',ms)
        end
    end
end

save([fnout 'eyeStruct'],'msEye')