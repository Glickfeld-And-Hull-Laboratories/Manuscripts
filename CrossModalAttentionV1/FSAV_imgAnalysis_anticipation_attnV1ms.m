clear all
close all
ds = 'FSAV_attentionV1';
cellsOrDendrites = 1;
doLoadPreviousAnalysis = true;
%%
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms

eval(ds)
titleStr = ds(6:end);
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

if cellsOrDendrites == 1
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));
    fnout = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        [titleStr '_anticipation_']); 
elseif cellsOrDendrites == 2
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));
    fnout = fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        [titleStr '_den_anticipation_']); 
end

if strcmp(ds,'FSAV_attentionV1')
    load(fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs',...
        'V1_100ms_naive_anticipation_imgStats'))
    imgStats_naive = imgStats;
    clear imgStats
end

%%
nav = 2;
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nFrames1s = frameRateHz;
nexp = size(expt,2);
nCycles = 8;
% lateCycles = 5:nCycles;
lateWinFr = (45:88)+nBaselineFr;
firstWinFr = (3:44)+nBaselineFr;
minTargetRT = (nVisDelayFr_target+respwin_target(1)-nBaselineFr)./frameRateHz*1000;

oriBinSize = 45;
orientations = 0:oriBinSize:(180-oriBinSize);
oriBinEdges = [0, (oriBinSize/2):oriBinSize:(180-(oriBinSize/2)), 180];
nOri = length(orientations);

nMovWin = 15;
movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelFr_target = 30:(30+nMovWin-1);
% movWinLabelMs = 

% minCellN_SIFRmatch = 36;

trOutType = {'h';'m';'fa';'cr'};
trOutTypeName = {'H-All';'H-HT';'H-ET';'M-All';'M-HT';'M-ET';'FA';'CR'};
%% pool experiment data
if doLoadPreviousAnalysis
    load([fnout 'imgAnalysisData'])
else
    antiDataExpt = struct;
    oriTuningExpt = struct;
    targetDataExpt = struct;
    decodeDataExpt = struct;
    decodeDataExpt.av(visualTrials).name = 'Visual';
    decodeDataExpt.av(auditoryTrials).name = 'Auditory';
    nTargets = 2; %sum(unique(targetInd) > 1);
    nTrialsPerExpt = nan(1,nexp);
    rewExptInd = nan(1,nexp);
    for imouse = 1:size(mouse,2)
        for iexp = 1:size(mouse(imouse).expt,2)
            if imouse == 1 && iexp == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end

            d = mouse(imouse).expt(iexp);
            exptID = strcmp({expt.SubNum},mouse(imouse).mouse_name) & strcmp({expt.date},d.date);
            if strcmp(ds,'FSAV_attentionV1')
                if expt(exptID).catchRew == 1
                    rewExptInd(exptN) = true;
                else
                    rewExptInd(exptN) = false;
                end
            else
                rewExptInd(exptN) = false;
            end
            nTrialsPerExpt(exptN) = length(d.av(visualTrials).align(1).outcome) + ...
                length(d.av(auditoryTrials).align(1).outcome);

            cycLengthFr = d.info.cycTimeFrames;
            nCycLongTC = ceil(longTrialLengthFr./cycLengthFr);

            maxCycles = max(cat(2,d.av(visualTrials).align(alignStart).nCycles,...
                d.av(auditoryTrials).align(alignStart).nCycles));

            tc_AV = [];
    %         cycTC = cell(2,nCycles);
            cycTC = cell(2,maxCycles);
            longTC = cell(2,1);
            for iav = 1:2
                dd = d.av(iav);
                trOut = [];
                trStim = [];
                trResp = [];
                trTC = [];
                trResp_movWin = [];
                for ialign = 2:4
                    ddd = dd.align(ialign);
    %                 if iav 
    %                     fprintf('%s/%s auditory trials < 200 ms\n', ...
    %                         num2str(sum(ddd.reactTime < 200)),num2str(length(ddd.reactTime)))
    %                 end
                    tc = ddd.respTC;
                    if ialign == 4
                        tc_temp = circshift(tc,-1,1);
                        tc_temp(end,:,:) = nan;
                        trTC = cat(3,trTC,tc);
                        trResp = cat(2,trResp,squeeze(mean(tc(respwin_target,:,:),1)...
                            - mean(tc(basewin_0_target,:,:),1)));
                    else
                        trResp = cat(2,trResp,squeeze(mean(tc(respwin,:,:),1)...
                            - mean(tc(basewin_0,:,:),1)));
                        trTC = cat(3,trTC,circshift(tc,-1,1)- mean(tc(basewin_0_target,:,:),1));
                    end
                    nfr = length(respwin)-1;
                    nMovWin = 15;
                    trResp_movWin_temp = nan(size(tc,2),size(tc,3),nMovWin);
                    if ialign == 4
                        tempwinstart = 28;
                    else
                        tempwinstart = 29;
                    end
                    for iwin = 1:nMovWin
                        tempwin = tempwinstart+(iwin-1):(tempwinstart+(iwin-1)+nfr);
                        if ialign == 4
                            trResp_movWin_temp(:,:,iwin) = squeeze(mean(tc(tempwin,:,:),1)) - ...
                                squeeze(mean(tc(basewin_0_target,:,:),1));
                        else
                            trResp_movWin_temp(:,:,iwin) = squeeze(mean(tc(tempwin,:,:),1)) - ...
                                squeeze(mean(tc(basewin_0,:,:),1));
                        end
                    end
                    trResp_movWin = cat(2,trResp_movWin,trResp_movWin_temp);
                    if ialign == alignCR
                        trOut = cat(2,trOut,repmat({'cr'},[1,length(ddd.outcome)]));
                    else
                        trOut = cat(2,trOut,ddd.outcome);
                    end
                    if isempty(ddd.ori) & isempty(ddd.amp)
                        trStim = cat(2,trStim,zeros(1,length(ddd.outcome)));
                    elseif iav == 1
                        trStim = cat(2,trStim,ddd.ori);
                    elseif iav == 2
                        trStim = cat(2,trStim,ddd.amp);
                    end
                end
                trOut(strcmp(trOut,'success')) = {'h'};
                trOut(strcmp(trOut,'ignore')) = {'m'};
                trOut(strcmp(trOut,'failure')) = {'fa'};

                decodeDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
                decodeDataExpt(exptN).av(iav).outcome = trOut;
                decodeDataExpt(exptN).av(iav).stim = trStim;
                decodeDataExpt(exptN).av(iav).resp = trResp;
                decodeDataExpt(exptN).av(iav).movWinResp = trResp_movWin;
                decodeDataExpt(exptN).av(iav).tc = trTC(1:(nBaselineFr*2),:,:);

                if iav == 1
                    binEdges = oriBins;
                elseif iav == 2
                    binEdges = ampBins;
                end
                trStimID = discretize(trStim,binEdges);
                rng(0)
                matchTrInd = getMatchedOutcomeTrialIndex(trOut,trStimID,minTrN);
                trOutTC = cell(1,length(trOutTypeName));
                for iout = 1:length(trOutType)
                    if iout == 1 || iout == 2
                        for istim = 2:3
                            ind = trStimID == istim;
                            if istim == 2
                                ind2 = ((iout-1)*3)+2;
                            elseif istim == 3
                                ind2 = ((iout-1)*3)+3;
                            end
                            ind3 = ind & strcmp(trOut,trOutType(iout));
                            if sum(ind3) == 0
                                trOutTC{ind2} = nan(size(trTC,1),size(trTC,2),1);
                            else
                                trOutTC{ind2} = trTC(:,:,ind3);
                            end
    %                         fprintf(num2str(ind2))
                        end
                        ind2 = ((iout-1)*3)+1;
                        trOutTC{ind2} = trTC(:,:,strcmp(trOut,trOutType(iout)) & matchTrInd);
    %                     fprintf(num2str(ind2))
                    else
    %                     fprintf(num2str(iout+4))
                        trOutTC{iout+4} = trTC(:,:,strcmp(trOut,trOutType(iout)));
                    end
                end
                targetDataExpt(exptN).av(iav).trOutTC = trOutTC;


                de = d.av(iav).align(alignStart);
                if strcmp(ds,'FSAV_attentionV1')
                    hits = strcmp(de.outcome,'success');
                elseif strcmp(ds,'FSAV_V1_100ms_naive')
                    hits = true(1,length(de.outcome));
                end
                misses = strcmp(de.outcome,'ignore');
                tc_AV = cat(3,tc_AV,de.respTC(:,:,hits)); % make hits or misses
    %             for icyc = 1:nCycles
                for icyc = 1:maxCycles
                    tc = de.respTC(:,:,de.nCycles >= icyc & (hits | misses));
                    cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                    cycTC{iav,icyc} = tc(...
                        (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
                end
                longTC{iav} = de.respTC(1:(longTrialLengthFr+nBaselineFr),:,...
                    de.nCycles >= nCycLongTC & (hits | misses));

                de = d.av(iav).align(alignTarget);
                if iav == 1
                    binEdges = oriBins;
                    targets = de.ori;
                elseif iav == 2
                    binEdges = ampBins;
                    targets = de.amp;
                end
                tc = de.respTC(1:(nBaselineFr*2),:,:);
                targetInd = discretize(targets,binEdges);
                targetDataExpt(exptN).av(iav).tc = cell(1,2);
                for itar = 1:nTargets
                    ind = targetInd == itar+1;
                    targetDataExpt(exptN).av(iav).tc{itar} = tc(:,:,ind);
                    targetDataExpt(exptN).av(iav).targets{itar} = targets(ind);
                end
            end
            antiDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            antiDataExpt(exptN).exptCycLengthFr = cycLengthFr;
            antiDataExpt(exptN).longTC = longTC;
            antiDataExpt(exptN).cycTC = cycTC;
            targetDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            if size(d.av(visualTrials).align,2) > 4
                r = squeeze(mean(...
                    d.av(visualTrials).align(5).respTC(respwin_target,:,:),1) - ...
                    mean(d.av(visualTrials).align(5).respTC(basewin_0_target,:,:),1));
                trOut =  d.av(visualTrials).align(5).outcome;
                ind = cellfun(@(x) ~isempty(x),trOut);
                trOut = trOut(ind);
                r = r(:,ind);
                trOut(strcmp(trOut,'FA')) = {'h'};
                trOut(strcmp(trOut,'CR')) = {'m'};
                decodeDataExpt(exptN).av(visualTrials).catchResp = r;
                decodeDataExpt(exptN).av(visualTrials).catchOutcome = trOut;
                decodeDataExpt(exptN).av(visualTrials).catchStim = ...
                    d.av(visualTrials).align(5).ori;
            end
                   
            
            do = d.oriTuning;
            if size(do.oriResp,2) == 8
                oriTuningExpt(exptN).oriResp = do.oriResp(:,1:2:8);
                oriTuningExpt(exptN).oriRespErr = do.oriRespSem(:,1:2:8);
            elseif size(do.oriResp,2) ~= 4
                error('error in orientations used for passivie tuning')
            else
                oriTuningExpt(exptN).oriResp = do.oriResp;
                oriTuningExpt(exptN).oriRespErr = do.oriRespSem;
            end
            oriTuningExpt(exptN).fit = do.oriFit;
            oriTuningExpt(exptN).isTuned = do.oriFitReliability < tuningReliabilityThresh;
            [~,oriTuningExpt(exptN).fitPeak] = max(do.oriFit,[],1);
            [~,oriPref] = histc(oriTuningExpt(exptN).fitPeak,oriBinEdges);
            oriPref(oriPref == length(orientations)+1 | oriPref == length(oriBinEdges) == 1) = 1;
            oriTuningExpt(exptN).oriPref = oriPref;
            oriTuningExpt(exptN).tuningReliability = do.oriFitReliability;
        end
    end

    shortCycExptInd = cell2mat({antiDataExpt.exptCycLengthFr}) == 11;
    isShortCycExpt = [];
    respCellsExpt = struct;
    for iexp = 1:nexp
        firstTCAV = cat(3,antiDataExpt(iexp).cycTC{visualTrials,1},...
            antiDataExpt(iexp).cycTC{auditoryTrials,1});

        lateTC_vis = [];
        lateTC_aud = [];
        lateCycles = 5:size(antiDataExpt(iexp).cycTC,2);
        for icyc = 1:length(lateCycles)
            lateTC_vis = cat(3,lateTC_vis,...
                antiDataExpt(iexp).cycTC{visualTrials,lateCycles(icyc)});
            lateTC_aud = cat(3,lateTC_aud,...
                antiDataExpt(iexp).cycTC{auditoryTrials,lateCycles(icyc)});
        end
        lateTCAV = cat(3,lateTC_vis,lateTC_aud);

        longTCAV = cat(3,antiDataExpt(iexp).longTC{visualTrials,1},...
            antiDataExpt(iexp).longTC{auditoryTrials,1});

        firstRespCells = ttest(...
            squeeze(mean(firstTCAV(respwin,:,:),1)),...
            squeeze(mean(firstTCAV(basewin_0,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);
        lateRespCells = ttest(...
            squeeze(mean(longTCAV(lateWinFr,:,:),1)),...
            squeeze(mean(longTCAV(basewin,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);
        lateSuppCells = ttest(...
            squeeze(mean(longTCAV(lateWinFr,:,:),1)),...
            squeeze(mean(longTCAV(basewin,:,:),1)),...
            'dim',2,'tail','left','alpha',cellGroupsAlpha);
        lateCycRespCells = ttest(...
            squeeze(mean(lateTCAV(respwin,:,:),1)),...
            squeeze(mean(lateTCAV(basewin_0,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);

        lateCycRespCutoffPass = mean(mean(lateTCAV(respwin,:,:),3),1) > minRespThreshold_decode;

        eaTarRespCells = sum(cell2mat(cellfun(@(x) ...
            ttest(squeeze(mean(x(respwin_target,:,:),1)),...
            squeeze(mean(x(basewin_0_target,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha),...
            targetDataExpt(iexp).av(visualTrials).tc,'unif',0)),2) > 0;


        allTargetTC = [];
        for itar = 1:length(targetDataExpt(iexp).av(visualTrials).tc)
            allTargetTC = cat(3,allTargetTC,...
                targetDataExpt(iexp).av(visualTrials).tc{itar});
        end
        allTarRespCells = ttest(squeeze(mean(allTargetTC(respwin_target,:,:),1)),...
            squeeze(mean(allTargetTC(basewin_0_target,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);

        eaTarRespCutoffPass = sum(cell2mat(cellfun(@(x) mean(mean(x(respwin_target,:,:),3),1) - mean(mean(x(basewin_0_target,:,:),3),1),...
            targetDataExpt(iexp).av(visualTrials).tc,'unif',0)') > minRespThreshold_decode,1) > 0 ...
            | (mean(mean(allTargetTC(respwin_target,:,:),3),1) - mean(mean(allTargetTC(basewin_0_target,:,:),3),1)) > minRespThreshold_decode; 

        respCellsExpt(iexp).exptName = antiDataExpt(iexp).exptName;
        respCellsExpt(iexp).firstRespCells = firstRespCells;
        respCellsExpt(iexp).lateRespCells = lateRespCells;
        respCellsExpt(iexp).lateSuppCells = lateSuppCells;
        respCellsExpt(iexp).lateCycRespCells = lateCycRespCells;
        respCellsExpt(iexp).targetRespCells = eaTarRespCells | allTarRespCells;
    %     respCellsExpt(iexp).decodeAnalysisCells = ...
    %         (firstRespCells & mean(mean(firstTCAV(respwin,:,:),3),1)' > minRespThreshold)...
    %         | ((eaTarRespCells | allTarRespCells) & eaTarRespCutoffPass');
        respCellsExpt(iexp).decodeAnalysisCells = ...
            (lateCycRespCells & lateCycRespCutoffPass')...
            | ((eaTarRespCells | allTarRespCells) & eaTarRespCutoffPass');
        if shortCycExptInd(iexp)
            isShortCycExpt = cat(1,isShortCycExpt,true(length(firstRespCells),1));
        else
            isShortCycExpt = cat(1,isShortCycExpt,false(length(firstRespCells),1));
        end

    end

    if strcmp(ds,'FSAV_attentionV1')
%         dcInfo = struct;
        decodeAnalysis = struct;
        decodeAnalysis.av(visualTrials).name = 'Visual';
        decodeAnalysis.av(auditoryTrials).name = 'Auditory';
        for iexp = 1:nexp
            rng(0) % cells randomized and trials randomized
            cellInd = respCellsExpt(iexp).decodeAnalysisCells & ...
                oriTuningExpt(iexp).tuningReliability' <= tuningReliabilityThresh_decode;
            decodeAnalysis(iexp).nCells = sum(cellInd);
            if sum(cellInd) > maxCellN
                ind = find(cellInd);
                cellSampleID = randsample(ind,maxCellN);
                cellInd = false(length(cellInd),1);
                cellInd(cellSampleID) = true;
            end
            decodeAnalysis(iexp).nCellsSelected = sum(cellInd);
            decodeAnalysis(iexp).cellInd = cellInd;
%             dcInfo(iexp).cellInd = cellInd;
            respOther = cell(1,2);
            trOutOther = cell(1,2);
            trSampleIndOther = cell(1,2);
            detectGLMOther = cell(1,2);
            targetGLMOther = cell(1,2);
            if isfield(decodeDataExpt(iexp).av(visualTrials),'catchOutcome')
                respAllCells_av_withInv = zscore(cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).resp,...
                    decodeDataExpt(iexp).av(auditoryTrials).resp,...
                    decodeDataExpt(iexp).av(visualTrials).catchResp)');
                nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
                naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
                ninv = size(decodeDataExpt(iexp).av(visualTrials).catchResp,2);
                respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
                respAllCells_aud = respAllCells_av_withInv((nvis+1):(nvis+naud),:);
                catchResp = respAllCells_av_withInv((end-ninv+1):end,:);
            else
                respAllCells_av_withInv = zscore(cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).resp,...
                    decodeDataExpt(iexp).av(auditoryTrials).resp)');
                nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
                naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
                respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
                respAllCells_aud = respAllCells_av_withInv((nvis+1):end,:);
            end
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;
%                 respAllCells = zscore(decodeDataExpt(iexp).av(iav).resp');
                

                if iav == 1
                    respAllCells = respAllCells_vis;
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                elseif iav == 2
                    respAllCells = respAllCells_aud;
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                end
                nStimPerBin = histcounts(trStimID);
                minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
                respStimSort = cell(1,nStimBins);
                trOutStimSort = cell(1,nStimBins);
                for istim = 1:nStimBins
                    ind = find(trStimID == istim);
                    if length(ind) >= minTrN_mdl
                        if istim == 1
                            matchTrialsInd = [];
                            if sum(nStimPerBin >= minTrN_mdl) == 2
                                n = minBinN;
                            elseif minBinN == nStimPerBin(istim)
                                error('not enough FA/CR trials')
                            else
                                n = (nStimBins-1).*minBinN;
                                if n > length(ind)
                                    error('not enough FA/CR trials')
                                end
                            end
                            indSample = randsample(ind,n);
                            matchTrialsInd = cat(2,matchTrialsInd,indSample);
                        else
                            indSample = randsample(ind,minBinN);
                            matchTrialsInd = cat(2,matchTrialsInd,indSample);
                        end
                        respStimSort{istim} = respAllCells(indSample,cellInd);
                        trOutStimSort{istim} = trOut(indSample);
                    end
                end
                nMatchedTrials = cumsum(cellfun(@length,trOutStimSort));
                for istim = 1:nStimBins
                    if istim == 1
                        stimSortInd = cell(1,nStimBins);
                        stimSortInd{istim} = 1:nMatchedTrials;
                    else
                        stimSortInd{istim} = ...
                            (nMatchedTrials(istim-1)+1):nMatchedTrials(istim);
                    end
                end
                trSampleIndOther{iav} = matchTrialsInd;
%                 dcInfo(iexp).av(iav).resp = respAllCells(matchTrialsInd,:);
%                 dcInfo(iexp).av(iav).trOut = trOut(matchTrialsInd);
        %         matchTrialsInd = 1:size(respAllCells,2);
        %         resp = zscore(respAllCells(cellInd,matchTrialsInd)');
                decodeAnalysis(iexp).av(iav).respAllCells = respAllCells(matchTrialsInd,:);
                decodeAnalysis(iexp).av(iav).trOut = trOut(matchTrialsInd);
                resp = respAllCells(matchTrialsInd,cellInd);
        %         [c,score,latent] = pca(decodeDataExpt(iexp).av(iav).resp(:,matchTrialsInd)');
        %         resp_pcs = zscore(score(:,1:10));
        %         [detectTrIndAll, targetTrIndAll] = getStimAndBehaviorYs(trOut);
                [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));

                detectCorr = corr(detectTrInd,resp);
                targetCorr = corr(targetTrInd,resp);

                C = eye(size(resp,2));
                p=1;
                [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
                [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');

                fprintf('Expt %s, hold-out analysis\n',num2str(iexp))
                detectWeight = detectGLM.beta(2:end);
                targetWeight = targetGLM.beta(2:end);

                dv_detect = mean(detectTrInd);
                dv_target = mean(targetTrInd);
                
                pctCorrectDetect_train = getPctCorr_trainData(detectGLM,resp,detectTrInd,dv_detect);
                pctCorrectDetect_ho = getPctCorr_hoData(resp,detectTrInd,dv_detect);

                pctCorrectTarget_train = getPctCorr_trainData(targetGLM,resp,targetTrInd,dv_target);
                pctCorrectTarget_ho = getPctCorr_hoData(resp,targetTrInd,dv_target);
 
                pctCorrDetect_xStim_train = nan(1,nStimBins);
                pctCorrDetect_xStim_ho = nan(1,nStimBins);
                pctCorrTarget_xStim_train = nan(1,nStimBins);
                pctCorrTarget_xStim_ho = nan(1,nStimBins);
                for istim = 1:nStimBins
                    if nStimPerBin(istim) >= minTrN_mdl
                        [detectStimInd, targetStimInd] = getStimAndBehaviorYs(...
                            trOutStimSort{istim});
                        pctCorrDetect_xStim_train(istim) = getPctCorr_trainData(...
                            detectGLM,respStimSort{istim},detectStimInd,dv_detect);
                        pctCorrDetect_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
                            resp,detectTrInd,stimSortInd{istim},dv_detect);
                        pctCorrTarget_xStim_train(istim) = getPctCorr_trainData(...
                            targetGLM,respStimSort{istim},targetStimInd,dv_target);
                        pctCorrTarget_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
                            resp,targetTrInd,stimSortInd{istim},dv_target);
                    end
                end

                fprintf('Expt %s, starting resp win analysis\n',num2str(iexp))
                nwins = size(decodeDataExpt(iexp).av(iav).movWinResp,3);
                pctCorrectDetect_train_respwin = nan(1,nwins);
                pctCorrectDetect_ho_respwin = nan(1,nwins);
                pctCorrectTarget_train_respwin = nan(1,nwins);
                pctCorrectTarget_ho_respwin = nan(1,nwins);
                for iwin = 1:nwins
                    rAllCells = zscore(decodeDataExpt(iexp).av(iav).movWinResp(:,:,iwin)');
                    r = rAllCells(matchTrialsInd,cellInd);
                    [~,~,detectGLM_temp] = glmfit(r*C,detectTrInd,'binomial');
                    [~,~,targetGLM_temp] = glmfit(r*C,targetTrInd,'binomial');

                    pctCorrectDetect_train_respwin(iwin) = getPctCorr_trainData(...
                        detectGLM_temp,r,detectTrInd,dv_detect);
                    pctCorrectDetect_ho_respwin(iwin) = getPctCorr_hoData(r,detectTrInd,dv_detect);
                    pctCorrectTarget_train_respwin(iwin) = getPctCorr_trainData(...
                        targetGLM_temp,r,targetTrInd,dv_target);
                    pctCorrectTarget_ho_respwin(iwin) = getPctCorr_hoData(r,targetTrInd,dv_target);
                end

                respOther{iav} = resp;
                trOutOther{iav} = trOut(matchTrialsInd);
                detectGLMOther{iav} = detectGLM;
                targetGLMOther{iav} = targetGLM;        

                decodeAnalysis(iexp).av(iav).dvDetect = dv_detect;
                decodeAnalysis(iexp).av(iav).dvTarget = dv_target;
                decodeAnalysis(iexp).av(iav).correlationDetect = detectCorr;
                decodeAnalysis(iexp).av(iav).correlationTarget = targetCorr;
                decodeAnalysis(iexp).av(iav).weightDetect = detectWeight;
                decodeAnalysis(iexp).av(iav).weightTarget = targetWeight;
                decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_train = pctCorrectDetect_train;
                decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_holdout = pctCorrectDetect_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_train = pctCorrectTarget_train;
                decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_holdout = pctCorrectTarget_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimDetect_train = pctCorrDetect_xStim_train;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimDetect_holdout = pctCorrDetect_xStim_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_train = pctCorrTarget_xStim_train;
                decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_holdout = pctCorrTarget_xStim_ho;
                decodeAnalysis(iexp).av(iav).pctCorrectDetectMovRespWin_train = pctCorrectDetect_train_respwin;
                decodeAnalysis(iexp).av(iav).pctCorrectDetectMovRespWin_holdout = pctCorrectDetect_ho_respwin;
                decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_train = pctCorrectTarget_train_respwin;
                decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_holdout = pctCorrectTarget_ho_respwin;
%                 decodeAnalysis(iexp).av(iav).validMatchedPctCorrectDetect_holdout = validCatchMatchDetect;
%                 decodeAnalysis(iexp).av(iav).validMatchedPctCorrectTarget_holdout = validCatchMatchTarget;
            end
            fprintf('Expt %s: Testing opposite model...\n',num2str(iexp))
            for iav = 1:2
                if iav == 1
                    otherAV = 2;
                    if isfield(decodeDataExpt(iexp).av(iav),'catchOutcome')
                        [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{otherAV});
                        [catchDetectInd,catchTargetInd] = getStimAndBehaviorYs(...
                            decodeDataExpt(iexp).av(iav).catchOutcome);
                        catchTargetResp = catchResp(:,cellInd);
                        nCatch = length(catchDetectInd);
                        catchDistRespAll = respOther{otherAV}(targetTrInd==0,:);
                        catchDistInd = randsample(sum(targetTrInd==0),nCatch);
                        catchRespBalanced = cat(1,catchDistRespAll(catchDistInd,:),...
                            catchTargetResp);
                        catchDetectIndBalanced = cat(1,detectTrInd(catchDistInd),...
                            catchDetectInd);
                        catchTargetIndBalanced = cat(1,targetTrInd(catchDistInd),...
                            catchTargetInd);

                        catchPctCorrectDetect = getPctCorr_trainData(...
                            detectGLMOther{iav},catchRespBalanced,catchDetectIndBalanced,...
                            decodeAnalysis(iexp).av(iav).dvDetect);
                        catchPctCorrectTarget = getPctCorr_trainData(...
                            targetGLMOther{iav},catchRespBalanced,catchTargetIndBalanced,...
                            decodeAnalysis(iexp).av(iav).dvTarget);
                        catchPctCorrectDetect_audModel = getPctCorr_trainData(...
                            detectGLMOther{auditoryTrials},catchRespBalanced,...
                            catchDetectIndBalanced,...
                            decodeAnalysis(iexp).av(auditoryTrials).dvDetect);            
                        catchPctCorrectTarget_audModel = getPctCorr_trainData(...
                            targetGLMOther{auditoryTrials},catchRespBalanced,...
                            catchTargetIndBalanced,...
                            decodeAnalysis(iexp).av(auditoryTrials).dvTarget);

                        trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                        matchedTrialsID = trStimID(trSampleIndOther{iav});
                        catchStimID = cat(2,ones(1,nCatch),...
                            discretize(decodeDataExpt(iexp).av(iav).catchStim,oriBins));

                        [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{iav});
                        catchMatchTrials = cell2mat(getMatchedValidTrialIndex(...
                            matchedTrialsID,catchStimID)); 
                        validCatchMatchDetect = getPctCorr_hoData_subGroup(...
                            respOther{iav},detectTrInd,catchMatchTrials,decodeAnalysis(iexp).av(iav).dvDetect);  
                        validCatchMatchTarget = getPctCorr_hoData_subGroup(...
                            respOther{iav},targetTrInd,catchMatchTrials,decodeAnalysis(iexp).av(iav).dvTarget);
                        
                        trStimID = discretize(decodeDataExpt(iexp).av(auditoryTrials).stim,ampBins);
                        matchedTrialsID = trStimID(trSampleIndOther{auditoryTrials});
                        matchedTrialsID(matchedTrialsID > 1) = 2;
                        catchStimID(catchStimID > 1) = 2;
                        [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{auditoryTrials});
                        catchMatchTrials = cell2mat(getMatchedValidTrialIndex(...
                            matchedTrialsID,catchStimID));
                        validCatchMatchDetect_aud = getPctCorr_hoData_subGroup(...
                            respOther{auditoryTrials},detectTrInd,catchMatchTrials,...
                            decodeAnalysis(iexp).av(auditoryTrials).dvDetect); 
                        validCatchMatchTarget_aud = getPctCorr_hoData_subGroup(...
                            respOther{auditoryTrials},targetTrInd,catchMatchTrials,...
                            decodeAnalysis(iexp).av(auditoryTrials).dvTarget);  
                       
                    else
                        validCatchMatchDetect = [];
                        validCatchMatchTarget = [];
                        validCatchMatchDetect_aud = [];
                        validCatchMatchTarget_aud = [];
                        catchPctCorrectDetect = [];
                        catchPctCorrectTarget = [];
                        catchPctCorrectDetect_audModel = [];
                        catchPctCorrectTarget_audModel = [];
                        
                    end
                elseif iav == 2
                    otherAV = 1;
                    validCatchMatchDetect = [];
                    validCatchMatchTarget = [];
                    validCatchMatchDetect_aud = [];
                    validCatchMatchTarget_aud = [];
                    catchPctCorrectDetect = [];
                    catchPctCorrectTarget = [];
                    catchPctCorrectDetect_audModel = [];
                    catchPctCorrectTarget_audModel = [];
                end
                resp = respOther{otherAV};
                [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{otherAV});

                dv_detect = mean(detectTrInd);
                dv_target = mean(targetTrInd);

                pctCorrectDetect = getPctCorr_trainData(detectGLMOther{iav},resp,detectTrInd,dv_detect);
                pctCorrectTarget = getPctCorr_trainData(targetGLMOther{iav},resp,targetTrInd,dv_target);
                
                stimIndTarget = targetTrInd == 0;
                pctCorrectDetectxStim = nan(1,2);
                pctCorrectDetectxStim(1) = getPctCorr_trainData(...
                    detectGLMOther{iav},resp(stimIndTarget,:),...
                    detectTrInd(stimIndTarget),dv_detect);
                pctCorrectDetectxStim(2) = getPctCorr_trainData(...
                    detectGLMOther{iav},resp(~stimIndTarget,:),...
                    detectTrInd(~stimIndTarget),dv_detect);
                pctCorrectTargetxStim = nan(1,2);
                pctCorrectTargetxStim(1) = getPctCorr_trainData(...
                    targetGLMOther{iav},resp(stimIndTarget,:),...
                    targetTrInd(stimIndTarget),dv_target);
                pctCorrectTargetxStim(2) = getPctCorr_trainData(...
                    targetGLMOther{iav},resp(~stimIndTarget,:),...
                    targetTrInd(~stimIndTarget),dv_target);
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_otherAV = pctCorrectDetect;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_otherAV = pctCorrectTarget;
                decodeAnalysis(iexp).av(iav).pctCorrectDetectxStim_otherAV = pctCorrectDetectxStim;
                decodeAnalysis(iexp).av(iav).pctCorrectTargetxStim_otherAV = pctCorrectTargetxStim;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectDetect_holdout = validCatchMatchDetect;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectTarget_holdout = validCatchMatchTarget;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectDetect = catchPctCorrectDetect;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectTarget = catchPctCorrectTarget;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectDetect_testAudModel = catchPctCorrectDetect_audModel;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectTarget_testAudModel = catchPctCorrectTarget_audModel;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectDetect_testAudModel = validCatchMatchDetect_aud;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectTarget_testAudModel = validCatchMatchTarget_aud;
            end
            rng(0)
            randInd = cellfun(@(x) randsample(length(x),floor(length(x)/2)),...
                trOutOther,'unif',0);
            otherRandInd = cellfun(@(x,y) setdiff(1:length(x),y),trOutOther,randInd,'unif',0);
            nt = cellfun(@length,randInd);
            nt_other = cellfun(@length,otherRandInd);
%             nt = cellfun(@length,trOutOther);
            respAV = cat(1,respOther{visualTrials}(randInd{visualTrials},:),...
                respOther{auditoryTrials}(randInd{auditoryTrials},:));
%             respAV = cat(1,respOther{visualTrials},...
%                 respOther{auditoryTrials});
            [detectTrInd, targetTrInd] = getStimAndBehaviorYs(cat(2,...
                trOutOther{visualTrials}(randInd{visualTrials}),trOutOther{auditoryTrials}(randInd{auditoryTrials})));
            [detectTrInd_other, targetTrInd_other] = getStimAndBehaviorYs(cat(2,...
                trOutOther{visualTrials}(otherRandInd{visualTrials}),trOutOther{auditoryTrials}(otherRandInd{auditoryTrials})));
%             [detectTrInd, targetTrInd] = getStimAndBehaviorYs(cat(2,...
%                 trOutOther{visualTrials},trOutOther{auditoryTrials}));
            dv_detect = mean(detectTrInd);
            dv_target = mean(targetTrInd);
            C = eye(size(respAV,2));
            p=1;
            [~,~,detectGLM] = glmfit(respAV*C,detectTrInd,'binomial');
            [~,~,targetGLM] = glmfit(respAV*C,targetTrInd,'binomial');
            pctCorrectDetect_combo_vis = cmbPctCorrect(nt(visualTrials),...
                getPctCorr_hoData_subGroup(respAV,detectTrInd,1:nt(visualTrials),dv_detect),...
                nt_other(visualTrials),...
                getPctCorr_trainData(detectGLM,respOther{visualTrials}(otherRandInd{visualTrials},:),...
                    detectTrInd_other(1:nt_other),dv_detect));
            pctCorrectTarget_combo_vis = cmbPctCorrect(nt(visualTrials),...
                getPctCorr_hoData_subGroup(respAV,targetTrInd,1:nt(visualTrials),dv_target),...
                nt_other(visualTrials),...
                getPctCorr_trainData(targetGLM,respOther{visualTrials}(otherRandInd{visualTrials},:),...
                    targetTrInd_other(1:nt_other),dv_target));
            pctCorrectDetect_combo_aud = cmbPctCorrect(nt(auditoryTrials),...
                getPctCorr_hoData_subGroup(respAV,detectTrInd,(nt(visualTrials)+1):sum(nt),dv_detect),...
                nt_other(auditoryTrials),...
                getPctCorr_trainData(detectGLM,respOther{auditoryTrials}(otherRandInd{auditoryTrials},:),...
                    detectTrInd_other((nt(visualTrials)+1):sum(nt)),dv_detect));
            pctCorrectTarget_combo_aud = cmbPctCorrect(nt(auditoryTrials),...
                getPctCorr_hoData_subGroup(respAV,targetTrInd,(nt(visualTrials)+1):sum(nt),dv_target),...
                nt_other(auditoryTrials),...
                getPctCorr_trainData(targetGLM,respOther{auditoryTrials}(otherRandInd{auditoryTrials},:),...
                    targetTrInd_other((nt(visualTrials)+1):sum(nt)),dv_target));
                
%             pctCorrectDetect_combo_vis = getPctCorr_hoData_subGroup(respAV,detectTrInd,1:nt(visualTrials),dv_detect);
%             pctCorrectTarget_combo_vis = getPctCorr_hoData_subGroup(respAV,targetTrInd,1:nt(visualTrials),dv_target);
%             pctCorrectDetect_combo_aud = getPctCorr_hoData_subGroup(respAV,detectTrInd,(nt(visualTrials)+1):sum(nt),dv_detect);
%             pctCorrectTarget_combo_aud = getPctCorr_hoData_subGroup(respAV,targetTrInd,(nt(visualTrials)+1):sum(nt),dv_target);
            
            [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{visualTrials});
            pctCorrectDetect_vis = getPctCorr_hoData_subGroup(...
                respOther{visualTrials},...
                detectTrInd,randInd{visualTrials},mean(detectTrInd));
            pctCorrectTarget_vis = getPctCorr_hoData_subGroup(...
                respOther{visualTrials},...
                targetTrInd,randInd{visualTrials},mean(targetTrInd));
            [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{auditoryTrials});
            pctCorrectDetect_aud = getPctCorr_hoData_subGroup(...
                respOther{auditoryTrials},...
                detectTrInd,randInd{auditoryTrials},mean(detectTrInd));
            pctCorrectTarget_aud = getPctCorr_hoData_subGroup(...
                respOther{auditoryTrials},...
                targetTrInd,randInd{auditoryTrials},mean(targetTrInd));
            
            
            decodeAnalysis(iexp).comboTrainWeightDetect = detectGLM.beta(2:end);
            decodeAnalysis(iexp).comboTrainWeightTarget = targetGLM.beta(2:end);
            decodeAnalysis(iexp).av(visualTrials).pctCorrectDetect_comboTrain = pctCorrectDetect_combo_vis;
            decodeAnalysis(iexp).av(visualTrials).pctCorrectTarget_comboTrain = pctCorrectTarget_combo_vis;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectDetect_comboTrain = pctCorrectDetect_combo_aud;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectTarget_comboTrain = pctCorrectTarget_combo_aud;
            decodeAnalysis(iexp).av(visualTrials).pctCorrectDetect_visTrainComboMatch = pctCorrectDetect_vis;
            decodeAnalysis(iexp).av(visualTrials).pctCorrectTarget_visTrainComboMatch = pctCorrectTarget_vis;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectDetect_audTrainComboMatch = pctCorrectDetect_aud;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectTarget_audTrainComboMatch = pctCorrectTarget_aud;
        end

        nBoot = 1000;
        dcBootstrap = struct;
        mdlFail = cell(2,nexp);
        w_detect_combo_boot = cell(1,nexp);
        w_target_combo_boot = cell(1,nexp);
        for iexp = 1:nexp
            fprintf('Expt %s\n',num2str(iexp))
            rng(0)
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;
                respAllCells = zscore(decodeDataExpt(iexp).av(iav).resp');

                if iav == 1
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                elseif iav == 2
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                end
                nStimPerBin = histcounts(trStimID);
                minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
                for istim = 1:nStimBins
                    ind = find(trStimID == istim);
                    if length(ind) >= minTrN_mdl
                        if istim == 1
                            matchTrialsInd = [];
                            if sum(nStimPerBin >= minTrN_mdl) == 2
                                n = minBinN;
                            elseif minBinN == nStimPerBin(istim)
                                error('not enough FA/CR trials')
                            else
                                n = (nStimBins-1).*minBinN;
                                if n > length(ind)
                                    error('not enough FA/CR trials')
                                end
                            end
                            indSample = randsample(ind,n);
                            matchTrialsInd = cat(2,matchTrialsInd,indSample);
                        else
                            indSample = randsample(ind,minBinN);
                            matchTrialsInd = cat(2,matchTrialsInd,indSample);
                        end
                    end
                end
                mdlFail{iav,iexp} = false(1,nBoot);
                for iboot = 1:nBoot
                    fprintf('.')

                    nTotalCells = length(respCellsExpt(iexp).decodeAnalysisCells);
                    cellInd = randsample(nTotalCells,maxCellN);

                    resp = respAllCells(matchTrialsInd,cellInd);
                    [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));

                    detectCorr = corr(detectTrInd,resp);
                    targetCorr = corr(targetTrInd,resp);

                    C = eye(size(resp,2));
                    p=1;
                    [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
                    [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');
                    glmWarning = lastwarn;
                    if ~isempty(glmWarning)
                        mdlFail{iav,iexp}(iboot) = true;
                    end
                    glmWarning = [];
                    lastwarn('')

                    detectWeight = detectGLM.beta(2:end);
                    targetWeight = targetGLM.beta(2:end);

                    dcBootstrap(iexp).av(iav).cellInd{iboot} = cellInd;
                    dcBootstrap(iexp).av(iav).targetCorrelation{iboot} = targetCorr;
                    dcBootstrap(iexp).av(iav).detectCorrelation{iboot} = detectCorr;
                    dcBootstrap(iexp).av(iav).targetWeight{iboot} = targetWeight;
                    dcBootstrap(iexp).av(iav).detectWeight{iboot} = detectWeight;
                end
            end
            dcBootstrap(iexp).nTotalCells = nTotalCells;
            rng(0)
            n = length(decodeAnalysis(iexp).cellInd);
            wd = nan(n,nBoot);
            wt = nan(n,nBoot);
            [detectTrInd, targetTrInd] = getStimAndBehaviorYs(...
                cat(2,decodeAnalysis(iexp).av(1).trOut,decodeAnalysis(iexp).av(2).trOut));
            r = cat(1,decodeAnalysis(iexp).av(1).resp,decodeAnalysis(iexp).av(2).resp);
            for iboot = 1:nBoot
                cellInd = randsample(n,maxCellN);
                resp = r(:,cellInd);

                C = eye(size(resp,2));
                p=1;
                [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
                [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');
                wd(cellInd,iboot) = detectGLM.beta(2:end);
                wt(cellInd,iboot) = targetGLM.beta(2:end);
            end
            w_detect_combo_boot{iexp} = wd;
            w_target_combo_boot{iexp} = wt;
        end
        
        cellfun(@(x) sum(~isnan(x),2),w_detect_combo_boot)
        
        targetWeightBoot = cell(1,2);
        detectWeightBoot = cell(1,2);
        nWeights_target = cell(nexp,2);
        nWeights_detect = cell(nexp,2);
        for iexp = 1:nexp
            n = dcBootstrap(iexp).nTotalCells;
            for iav = 1:2
                d = dcBootstrap(iexp).av(iav);
                bootMatrix_w_target = nan(nBoot,n);
                bootMatrix_w_detect = nan(nBoot,n);
                for iboot = 1:nBoot
                    if mdlFail{iav,iexp}(iboot)
                        continue
                    end
                        ind = d.cellInd{iboot};
                        bootMatrix_w_target(iboot,ind) = d.targetWeight{iboot};
                        bootMatrix_w_detect(iboot,ind) = d.detectWeight{iboot};
                end
                targetWeightBoot{iav} = cat(2,targetWeightBoot{iav},nanmean(bootMatrix_w_target));
                detectWeightBoot{iav} = cat(2,detectWeightBoot{iav},nanmean(bootMatrix_w_detect));
                
                nWeights_target{iexp,iav} = sum(~isnan(bootMatrix_w_target));
                nWeights_detect{iexp,iav} = sum(~isnan(bootMatrix_w_detect));
            end
        end
        
    nCellsSelected = cell2mat({decodeAnalysis.nCellsSelected});    
    fprintf('%s +/- %s (%s-%s) cells selected\n',num2str(mean(nCellsSelected)),...
        num2str(ste(nCellsSelected,2)),num2str(min(nCellsSelected)),num2str(max(nCellsSelected)))
    fprintf('%s +/- %s (%s-%s) target weights from boot\n', ...
        num2str(mean(cat(2,cell2mat(nWeights_target(:,1)'),cell2mat(nWeights_target(:,2)')))),...
        num2str(ste(cat(2,cell2mat(nWeights_target(:,1)'),cell2mat(nWeights_target(:,2)')),2)),...
        num2str(min(cat(2,cell2mat(nWeights_target(:,1)'),cell2mat(nWeights_target(:,2)')))),...
        num2str(max(cat(2,cell2mat(nWeights_target(:,1)'),cell2mat(nWeights_target(:,2)')))))
    fprintf('%s +/- %s (%s-%s) detect weights from boot\n', ...
        num2str(mean(cat(2,cell2mat(nWeights_detect(:,1)'),cell2mat(nWeights_detect(:,2)')))),...
        num2str(ste(cat(2,cell2mat(nWeights_detect(:,1)'),cell2mat(nWeights_detect(:,2)')),2)),...
        num2str(min(cat(2,cell2mat(nWeights_detect(:,1)'),cell2mat(nWeights_detect(:,2)')))),...
        num2str(max(cat(2,cell2mat(nWeights_detect(:,1)'),cell2mat(nWeights_detect(:,2)')))))
    else
        decodeAnalysis = struct;
        decodeAnalysis.av(visualTrials).name = 'Visual';
        decodeAnalysis.av(auditoryTrials).name = 'Auditory';
        for iexp = 1:nexp
            rng(0) % cells randomized and trials randomized
            cellInd = respCellsExpt(iexp).decodeAnalysisCells & ...
                oriTuningExpt(iexp).tuningReliability' <= tuningReliabilityThresh_decode;
            decodeAnalysis(iexp).nCells = sum(cellInd);
            if sum(cellInd) > maxCellN
                ind = find(cellInd);
                cellSampleID = randsample(ind,maxCellN);
                cellInd = false(length(cellInd),1);
                cellInd(cellSampleID) = true;
            elseif sum(cellInd) < minCellN
                 cellInd = respCellsExpt(iexp).decodeAnalysisCells;
                 if sum(cellInd) > maxCellN
                    ind = find(cellInd);
                    cellSampleID = randsample(ind,maxCellN);
                    cellInd = false(length(cellInd),1);
                    cellInd(cellSampleID) = true;
                 end
            end
            decodeAnalysis(iexp).nCellsSelected = sum(cellInd);
            decodeAnalysis(iexp).cellInd = cellInd;
            respOther = cell(1,2);
            trOutOther = cell(1,2);
    %         detectGLMOther = cell(1,2);
            targetGLMOther = cell(1,2);
            respAllCells_av_withInv = zscore(cat(2,...
                decodeDataExpt(iexp).av(visualTrials).resp,...
                decodeDataExpt(iexp).av(auditoryTrials).resp)');
            nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
            naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
            respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
            respAllCells_aud = respAllCells_av_withInv((nvis+1):end,:);
            for iav = 1:2
                if sum(cellInd) < minCellN
                    decodeAnalysis(iexp).av(iav).dvDetect = nan;
                    decodeAnalysis(iexp).av(iav).dvTarget = nan;
                    decodeAnalysis(iexp).av(iav).correlationDetect = [];
                    decodeAnalysis(iexp).av(iav).correlationTarget = [];
                    decodeAnalysis(iexp).av(iav).weightDetect = [];
                    decodeAnalysis(iexp).av(iav).weightTarget = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_train = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_holdout = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_train = nan;
                    decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_holdout = nan;
                    decodeAnalysis(iexp).av(iav).pctCorrectXStimDetect_train = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectXStimDetect_holdout = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_train = nan(1,3);
                    decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_holdout = nan(1,3);
                    decodeAnalysis(iexp).av(iav).pctCorrectDetectMovRespWin_train = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectDetectMovRespWin_holdout = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_train = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_holdout = [];

                else
                    trOut = decodeDataExpt(iexp).av(iav).outcome;
%                     respAllCells = zscore(decodeDataExpt(iexp).av(iav).resp');

                    if iav == 1
                        respAllCells = respAllCells_vis;
                        trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                    elseif iav == 2
                        respAllCells = respAllCells_aud;
                        trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                    end
                    trOut(trStimID == 1) = {'cr'};
                    trOut(trStimID ~= 1) = {'h'};
                    nStimPerBin = histcounts(trStimID);
                    minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
                    respStimSort = cell(1,nStimBins);
                    trOutStimSort = cell(1,nStimBins);
                    for istim = 1:nStimBins
                        ind = find(trStimID == istim);
                        if length(ind) >= minTrN_mdl
                            if istim == 1
                                matchTrialsInd = [];
                                if sum(nStimPerBin >= minTrN_mdl) == 2
                                    n = minBinN;
                                elseif minBinN == nStimPerBin(istim)
                                    error('not enough FA/CR trials')
                                else
                                    n = (nStimBins-1).*minBinN;
                                    if n > length(ind)
                                        error('not enough FA/CR trials')
                                    end
                                end
                                indSample = randsample(ind,n);
                                matchTrialsInd = cat(2,matchTrialsInd,indSample);
                            else
                                indSample = randsample(ind,minBinN);
                                matchTrialsInd = cat(2,matchTrialsInd,indSample);
                            end
                            respStimSort{istim} = respAllCells(indSample,cellInd);
                            trOutStimSort{istim} = trOut(indSample);
                        end
                    end
                    nMatchedTrials = cumsum(cellfun(@length,trOutStimSort));
                    for istim = 1:nStimBins
                        if istim == 1
                            stimSortInd = cell(1,nStimBins);
                            stimSortInd{istim} = 1:nMatchedTrials;
                        else
                            stimSortInd{istim} = ...
                                (nMatchedTrials(istim-1)+1):nMatchedTrials(istim);
                        end
                    end

            %         matchTrialsInd = 1:size(respAllCells,2);
            %         resp = zscore(respAllCells(cellInd,matchTrialsInd)');
                    resp = respAllCells(matchTrialsInd,cellInd);
            %         [c,score,latent] = pca(decodeDataExpt(iexp).av(iav).resp(:,matchTrialsInd)');
            %         resp_pcs = zscore(score(:,1:10));
            %         [detectTrIndAll, targetTrIndAll] = getStimAndBehaviorYs(trOut);
                    [~, targetTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));

        %             detectCorr = corr(detectTrInd,resp);
                    targetCorr = corr(targetTrInd,resp);

                    C = eye(size(resp,2));
                    p=1;
        %             [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
                    [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');

                    fprintf('Expt %s, hold-out analysis\n',num2str(iexp))
        %             detectWeight = detectGLM.beta(2:end);
                    targetWeight = targetGLM.beta(2:end);

        %             dv_detect = mean(detectTrInd);
                    dv_target = mean(targetTrInd);

        %             pctCorrectDetect_train = getPctCorr_trainData(detectGLM,resp,detectTrInd,dv_detect);
        %             pctCorrectDetect_ho = getPctCorr_hoData(resp,detectTrInd,dv_detect);

                    pctCorrectTarget_train = getPctCorr_trainData(targetGLM,resp,targetTrInd,dv_target);
                    pctCorrectTarget_ho = getPctCorr_hoData(resp,targetTrInd,dv_target);

        %             pctCorrDetect_xStim_train = nan(1,nStimBins);
        %             pctCorrDetect_xStim_ho = nan(1,nStimBins);
                    pctCorrTarget_xStim_train = nan(1,nStimBins);
                    pctCorrTarget_xStim_ho = nan(1,nStimBins);
                    for istim = 1:nStimBins
                        if nStimPerBin(istim) >= minTrN_mdl
                            [~, targetStimInd] = getStimAndBehaviorYs(...
                                trOutStimSort{istim});
        %                     pctCorrDetect_xStim_train(istim) = getPctCorr_trainData(...
        %                         detectGLM,respStimSort{istim},detectStimInd,dv_detect);
        %                     pctCorrDetect_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
        %                         resp,detectTrInd,stimSortInd{istim},dv_detect);
                            pctCorrTarget_xStim_train(istim) = getPctCorr_trainData(...
                                targetGLM,respStimSort{istim},targetStimInd,dv_target);
                            pctCorrTarget_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
                                resp,targetTrInd,stimSortInd{istim},dv_target);
                        end
                    end

        %             fprintf('Expt %s, starting resp win analysis\n',num2str(iexp))
        %             nwins = size(decodeDataExpt(iexp).av(iav).movWinResp,3);
        % %             pctCorrectDetect_train_respwin = nan(1,nwins);
        % %             pctCorrectDetect_ho_respwin = nan(1,nwins);
        %             pctCorrectTarget_train_respwin = nan(1,nwins);
        %             pctCorrectTarget_ho_respwin = nan(1,nwins);
        %             for iwin = 1:nwins
        %                 rAllCells = zscore(decodeDataExpt(iexp).av(iav).movWinResp(:,:,iwin)');
        %                 r = rAllCells(matchTrialsInd,cellInd);
        % %                 [~,~,detectGLM_temp] = glmfit(r*C,detectTrInd,'binomial');
        %                 [~,~,targetGLM_temp] = glmfit(r*C,targetTrInd,'binomial');
        % 
        % %                 pctCorrectDetect_train_respwin(iwin) = getPctCorr_trainData(...
        % %                     detectGLM_temp,r,detectTrInd,dv_detect);
        % %                 pctCorrectDetect_ho_respwin(iwin) = getPctCorr_hoData(r,detectTrInd,dv_detect);
        %                 pctCorrectTarget_train_respwin(iwin) = getPctCorr_trainData(...
        %                     targetGLM_temp,r,targetTrInd,dv_target);
        %                 pctCorrectTarget_ho_respwin(iwin) = getPctCorr_hoData(r,targetTrInd,dv_target);
        %             end

                    respOther{iav} = resp;
                    trOutOther{iav} = trOut(matchTrialsInd);
        %             detectGLMOther{iav} = detectGLM;
                    targetGLMOther{iav} = targetGLM;        

                    decodeAnalysis(iexp).av(iav).dvDetect = nan;
                    decodeAnalysis(iexp).av(iav).dvTarget = dv_target;
                    decodeAnalysis(iexp).av(iav).correlationDetect = [];
                    decodeAnalysis(iexp).av(iav).correlationTarget = targetCorr;
                    decodeAnalysis(iexp).av(iav).weightDetect = [];
                    decodeAnalysis(iexp).av(iav).weightTarget = targetWeight;
                    decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_train = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_holdout = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_train = pctCorrectTarget_train;
                    decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_holdout = pctCorrectTarget_ho;
                    decodeAnalysis(iexp).av(iav).pctCorrectXStimDetect_train = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectXStimDetect_holdout = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_train = pctCorrTarget_xStim_train;
                    decodeAnalysis(iexp).av(iav).pctCorrectXStimTarget_holdout = pctCorrTarget_xStim_ho;
                    decodeAnalysis(iexp).av(iav).pctCorrectDetectMovRespWin_train = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectDetectMovRespWin_holdout = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_train = [];
                    decodeAnalysis(iexp).av(iav).pctCorrectTargetMovRespWin_holdout = [];
                end
            end
            fprintf('Expt %s: Testing opposite model...\n',num2str(iexp))
            for iav = 1:2
                if sum(cellInd) < minCellN
                    decodeAnalysis(iexp).av(iav).pctCorrectTarget_otherAV = nan;
                else
                    if iav == 1
                        otherAV = 2;
                    elseif iav == 2
                        otherAV = 1;
                    end
                    resp = respOther{otherAV};
                    [~, targetTrInd] = getStimAndBehaviorYs(trOutOther{otherAV});

        %             dv_detect = mean(detectTrInd);
                    dv_target = mean(targetTrInd);

        %             pctCorrectDetect = getPctCorr_trainData(detectGLMOther{iav},resp,detectTrInd,dv_detect);
                    pctCorrectTarget = getPctCorr_trainData(targetGLMOther{iav},resp,targetTrInd,dv_target);

        %             decodeAnalysis(iexp).av(iav).pctCorrectDetect_otherAV = pctCorrectDetect;
                    decodeAnalysis(iexp).av(iav).pctCorrectTarget_otherAV = pctCorrectTarget;
                end
            end
        end

        nBoot = 1000;
        dcBootstrap = struct;
        mdlFail = cell(2,nexp);
        w_detect_combo_boot = cell(1,nexp);
        w_target_combo_boot = cell(1,nexp);
        for iexp = 1:nexp
            fprintf('Expt %s\n',num2str(iexp))
            rng(0)
            nTotalCells = length(respCellsExpt(iexp).decodeAnalysisCells);
            for iav = 1:2
                if decodeAnalysis(iexp).nCellsSelected < minCellN
                    dcBootstrap(iexp).av(iav).targetWeight = cell(1,nBoot);
                else
                    trOut = decodeDataExpt(iexp).av(iav).outcome;
                    respAllCells = zscore(decodeDataExpt(iexp).av(iav).resp');

                    if iav == 1
                        trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                    elseif iav == 2
                        trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                    end
                    trOut(trStimID == 1) = {'cr'};
                    trOut(trStimID ~= 1) = {'h'};
                    nStimPerBin = histcounts(trStimID);
                    minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
                    for istim = 1:nStimBins
                        ind = find(trStimID == istim);
                        if length(ind) >= minTrN_mdl
                            if istim == 1
                                matchTrialsInd = [];
                                if sum(nStimPerBin >= minTrN_mdl) == 2
                                    n = minBinN;
                                elseif minBinN == nStimPerBin(istim)
                                    error('not enough FA/CR trials')
                                else
                                    n = (nStimBins-1).*minBinN;
                                    if n > length(ind)
                                        error('not enough FA/CR trials')
                                    end
                                end
                                indSample = randsample(ind,n);
                                matchTrialsInd = cat(2,matchTrialsInd,indSample);
                            else
                                indSample = randsample(ind,minBinN);
                                matchTrialsInd = cat(2,matchTrialsInd,indSample);
                            end
                        end
                    end
                    mdlFail{iav,iexp} = false(1,nBoot);
                    for iboot = 1:nBoot
                        fprintf('.')

                        nTotalCells = length(respCellsExpt(iexp).decodeAnalysisCells);
                        cellInd = randsample(nTotalCells,maxCellN);

                        resp = respAllCells(matchTrialsInd,cellInd);
                        [~, targetTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));

                        targetCorr = corr(targetTrInd,resp);

                        C = eye(size(resp,2));
                        p=1;
                        [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');
                        glmWarning = lastwarn;
                        if ~isempty(glmWarning)
                            mdlFail{iav,iexp}(iboot) = true;
                        end
                        glmWarning = [];
                        lastwarn('')

                        targetWeight = targetGLM.beta(2:end);

                        dcBootstrap(iexp).av(iav).cellInd{iboot} = cellInd;
                        dcBootstrap(iexp).av(iav).targetCorrelation{iboot} = targetCorr;
                        dcBootstrap(iexp).av(iav).detectCorrelation{iboot} = [];
                        dcBootstrap(iexp).av(iav).targetWeight{iboot} = targetWeight;
                        dcBootstrap(iexp).av(iav).detectWeight{iboot} = [];
                    end
                end
            end
            dcBootstrap(iexp).nTotalCells = nTotalCells;
        end
        
        
        
        targetWeightBoot = cell(1,2);
        for iexp = 1:nexp
            n = dcBootstrap(iexp).nTotalCells;
            for iav = 1:2
                if decodeAnalysis(iexp).nCellsSelected < minCellN
                    targetWeightBoot{iav} = cat(2,targetWeightBoot{iav},nan(1,dcBootstrap(iexp).nTotalCells));
                else
                    d = dcBootstrap(iexp).av(iav);
                    bootMatrix_w_target = nan(nBoot,n);
                    for iboot = 1:nBoot
                        if mdlFail{iav,iexp}(iboot)
                            continue
                        end
                            ind = d.cellInd{iboot};
                            bootMatrix_w_target(iboot,ind) = d.targetWeight{iboot};
                    end
                    targetWeightBoot{iav} = cat(2,targetWeightBoot{iav},nanmean(bootMatrix_w_target));
                end
            end
        end

    end

    antiAnalysis = struct;
    antiAnalysis.longTC = cell(1,3);
    antiAnalysis.longTCErr = cell(1,3);
    antiAnalysis.lateCycTC = cell(1,3);
    antiAnalysis.lateCycTCErr = cell(1,3);
    antiAnalysis.lateCycTC_shuff = cell(1,3);
    antiAnalysis.firstCycTC = cell(1,3);
    antiAnalysis.firstCycTCErr = cell(1,3);
    antiAnalysis.lateCycSI = [];
    antiAnalysis.lateCycSI_shuff = [];    
    antiAnalysis.lateCycAVauROC = []; 
    antiAnalysis.lateCycAVauROC_test = [];
    antiAnalysis.lateWinSI = [];
    antiAnalysis.lateCycAV95CITest = [];
    antiAnalysis.lateWinAV95CITest = [];
    antiAnalysis.lateCycAVShuffTest = [];
    antiAnalysis.lateWinAVShuffTest = [];
    antiAnalysis.lateCycFracSI = [];
    targetAnalysis = struct;
    targetAnalysis.tc = cell(3,2);
    targetAnalysis.targets = cell(1,2);
    targetAnalysis.targetsAuROC = [];
    targetAnalysis.targetsAuROCTest = [];
    targetAnalysis.trOutTC = cell(length(trOutTypeName),2);
    firstStimAuROC = [];
    lateStimAuROC = [];
    firstStimAuROCTest = [];
    lateStimAuROCTest = [];
    firstStimAuROC_aud = [];
    lateStimAuROC_aud = [];
    firstStimAuROCTest_aud = [];
    lateStimAuROCTest_aud = [];
    taskTuningTest = [];
    taskTuningPref = [];
    nModCellsPerExpt = nan(1,nexp);
    for iexp = 1:nexp
        longTC_vis = antiDataExpt(iexp).longTC{visualTrials};
        longTC_aud = antiDataExpt(iexp).longTC{auditoryTrials};

        nCells = size(longTC_vis,2);

        firstCycResp_vis = antiDataExpt(iexp).cycTC{visualTrials,1};
        firstCycResp_aud = antiDataExpt(iexp).cycTC{auditoryTrials,1};

        cycTC_vis = antiDataExpt(iexp).cycTC(visualTrials,:);
        cycTC_aud = antiDataExpt(iexp).cycTC(auditoryTrials,:);

        lateCycTC_vis = [];
        lateCycTC_aud = [];
        lateCycles = 5:size(antiDataExpt(iexp).cycTC,2);
        for icyc = 1:length(lateCycles)
            lateCycTC_vis = cat(3,lateCycTC_vis,...
                cycTC_vis{lateCycles(icyc)} - mean(cycTC_vis{lateCycles(icyc)}(basewin_0,:,:),1));
            lateCycTC_aud = cat(3,lateCycTC_aud,...
                cycTC_aud{lateCycles(icyc)} - mean(cycTC_aud{lateCycles(icyc)}(basewin_0,:,:),1));
        end
        
%         lateCycAVTest = ttest2(squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
%             squeeze(mean(lateCycTC_aud(respwin,:,:),1))','dim',1,'alpha',0.05);
%         lateWinAVTest = ttest2(squeeze(mean(longTC_vis(lateWinFr,:,:),1))',...
%             squeeze(mean(longTC_aud(lateWinFr,:,:),1))','dim',1,'alpha',0.05);
        rng(0)
        [lateCycSI95CITest,lateCycAVShuffTest] = testSelectivityIndex(...
            squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
            squeeze(mean(lateCycTC_aud(respwin,:,:),1))');
        [lateWinSI95CITest,lateWinAVShuffTest] = testSelectivityIndex(...
            squeeze(mean(longTC_vis(lateWinFr,:,:),1))',...
            squeeze(mean(longTC_aud(lateWinFr,:,:),1))');
        nModCellsPerExpt(iexp) = sum(lateCycSI95CITest & respCellsExpt(iexp).lateCycRespCells);
        
        rng(0)
        lateCycSI = getSelectivityIndex(squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
            squeeze(mean(lateCycTC_aud(respwin,:,:),1))');
        lateWinSI = getSelectivityIndex(squeeze(mean(longTC_vis(lateWinFr,:,:),1))',...
            squeeze(mean(longTC_aud(lateWinFr,:,:),1))');
        [lateCycAVauROC, lateCycAVauROC_test] = getAVauROC(squeeze(mean(lateCycTC_aud(respwin,:,:),1))',...
            squeeze(mean(lateCycTC_vis(respwin,:,:),1))');
        if size(lateCycTC_vis,3) ~= size(lateCycTC_aud,3)
            if size(lateCycTC_vis,3) > size(lateCycTC_aud,3)
                randInd = randsample(size(lateCycTC_vis,3),size(lateCycTC_aud,3));
                lateCycRespAll = cat(1,squeeze(mean(lateCycTC_vis(respwin,:,randInd),1))',...
                    squeeze(mean(lateCycTC_aud(respwin,:,:),1))');
            else
                randInd = randsample(size(lateCycTC_aud,3),size(lateCycTC_vis,3));
                lateCycRespAll = cat(1,squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
                    squeeze(mean(lateCycTC_aud(respwin,:,randInd),1))');
            end
        else
            lateCycRespAll = cat(1,squeeze(mean(lateCycTC_vis(respwin,:,:),1))',...
                squeeze(mean(lateCycTC_aud(respwin,:,randInd),1))');                
        end
        
        [shuffResp_vis, shuffResp_aud] = twoRandGroups(lateCycRespAll);
        lateCycSI_shuff = getSelectivityIndex(shuffResp_vis,shuffResp_aud);
        v=mean(mean(lateCycTC_vis(respwin,:,:),3),1);
        v(v<0) = 0;
        a=mean(mean(lateCycTC_aud(respwin,:,:),3),1);
        a(a<0) = 0;
        lateCycFracSI = (v-a)./v;
        lateCycFracSI_aud = (a-v)./a;
        lateCycFracSI(a > v) = -lateCycFracSI_aud(a > v);

        firstCycTC_vis = antiDataExpt(iexp).cycTC{visualTrials,1};
        firstCycTC_aud = antiDataExpt(iexp).cycTC{auditoryTrials,1};

        firstVisResp = squeeze(mean(firstCycTC_vis(respwin,:,:),1));
        lateVisResp = squeeze(mean(lateCycTC_vis(respwin,:,:),1));

        tarVisResp = cellfun(@(x) ...
            squeeze(mean(x(respwin_target,:,:),1) - mean(x(basewin_0_target,:,:),1)),...
            cat(2,targetDataExpt(iexp).av(visualTrials).tc,...
            {cat(3,targetDataExpt(iexp).av(visualTrials).tc{1},...
            targetDataExpt(iexp).av(visualTrials).tc{2})}),'unif',0);
        auroc_first = nan(nCells,3);
        auroc_first_test = nan(nCells,3);
        auroc_late = nan(nCells,3);
        auroc_late_test = nan(nCells,3);
        auroc_target = nan(nCells,1);
        auroc_target_test = nan(nCells,1);
        for icell = 1:nCells
            auroc_target(icell) = roc_gh(tarVisResp{1}(icell,:),tarVisResp{2}(icell,:));
            auroc_target_test(icell) = ranksum(...
                tarVisResp{1}(icell,:),tarVisResp{2}(icell,:)) < 0.05;

            auroc_first(icell,:) = cellfun(@(x) ...
                roc_gh(firstVisResp(icell,:),x(icell,:)),tarVisResp);
            auroc_first_test(icell,:) = cellfun(@(x) ...
                ranksum(firstVisResp(icell,:),x(icell,:)),tarVisResp);
            auroc_late(icell,:) = cellfun(@(x) ...
                roc_gh(lateVisResp(icell,:),x(icell,:)),tarVisResp);
            auroc_late_test(icell,:) = cellfun(@(x) ...
                ranksum(lateVisResp(icell,:),x(icell,:)),tarVisResp);
        end

        firstAudResp = squeeze(mean(firstCycTC_aud(respwin,:,:),1));
        lateAudResp = squeeze(mean(lateCycTC_aud(respwin,:,:),1));
        tarAudResp = cell2mat(cellfun(@(x) ...
            squeeze(mean(x(respwin_target,:,:),1) - mean(x(basewin_0_target,:,:),1)),...
            {cat(3,targetDataExpt(iexp).av(auditoryTrials).tc{1},...
            targetDataExpt(iexp).av(auditoryTrials).tc{2})},'unif',0));

        auroc_first_aud = nan(nCells,1);
        auroc_first_aud_test = nan(nCells,1);
        auroc_late_aud = nan(nCells,1);
        auroc_late_aud_test = nan(nCells,1);
        for icell = 1:nCells        
            auroc_first_aud(icell,:) =  ...
                roc_gh(firstAudResp(icell,:),tarAudResp(icell,:));
            auroc_first_aud_test(icell,:) =  ...
                ranksum(firstAudResp(icell,:),tarAudResp(icell,:));
            auroc_late_aud(icell,:) = ...
                roc_gh(firstAudResp(icell,:),tarAudResp(icell,:));
            auroc_late_aud_test(icell,:) =  ...
                ranksum(firstAudResp(icell,:),tarAudResp(icell,:));
        end

        stimResp4anova = cat(2,firstVisResp,tarVisResp{1},tarVisResp{2});
        stimID4anova = cat(2,ones(1,size(firstVisResp,2)),...
            ones(1,size(tarVisResp{1},2))*2,ones(1,size(tarVisResp{2},2))*3);
        taskStimAnova = nan(1,nCells);
        for icell = 1:nCells
            taskStimAnova(icell) = anova1(stimResp4anova(icell,:),stimID4anova,'off');
        end
        [~,taskTuningID] = max(cat(2,mean(firstVisResp,2),...
            mean(tarVisResp{1},2),mean(tarVisResp{2},2)),[],2);

        antiAnalysis.longTC{visualTrials} = cat(2,antiAnalysis.longTC{visualTrials},...
            mean(longTC_vis,3));
        antiAnalysis.longTC{auditoryTrials} = cat(2,antiAnalysis.longTC{auditoryTrials},...
            mean(longTC_aud,3));
        antiAnalysis.longTC{allTrialsInd} = cat(2,antiAnalysis.longTC{allTrialsInd},...
            mean(cat(3,longTC_vis,longTC_aud),3));
        antiAnalysis.longTCErr{allTrialsInd} = cat(2,antiAnalysis.longTCErr{allTrialsInd},...
            ste(cat(3,longTC_vis,longTC_aud),3));
        antiAnalysis.lateCycTC{visualTrials} = cat(2,antiAnalysis.lateCycTC{visualTrials},...
            mean(lateCycTC_vis,3));
        antiAnalysis.lateCycTC{auditoryTrials} = cat(2,antiAnalysis.lateCycTC{auditoryTrials},...
            mean(lateCycTC_aud,3));
        antiAnalysis.lateCycTC{allTrialsInd} = cat(2,antiAnalysis.lateCycTC{allTrialsInd},...
            mean(cat(3,lateCycTC_vis,lateCycTC_aud),3));
        antiAnalysis.lateCycTCErr{visualTrials} = cat(2,antiAnalysis.lateCycTCErr{visualTrials},...
            ste(lateCycTC_vis,3));
        antiAnalysis.lateCycTCErr{auditoryTrials} = cat(2,antiAnalysis.lateCycTCErr{auditoryTrials},...
            ste(lateCycTC_aud,3));
        antiAnalysis.lateCycTCErr{allTrialsInd} = cat(2,antiAnalysis.lateCycTCErr{allTrialsInd},...
            ste(cat(3,lateCycTC_vis,lateCycTC_aud),3));
        antiAnalysis.lateCycSI = cat(2,antiAnalysis.lateCycSI,lateCycSI);
        antiAnalysis.lateCycSI_shuff = cat(2,antiAnalysis.lateCycSI_shuff,lateCycSI_shuff);
        antiAnalysis.lateCycAVauROC = cat(2,antiAnalysis.lateCycAVauROC,lateCycAVauROC);
        antiAnalysis.lateCycAVauROC_test = cat(2,antiAnalysis.lateCycAVauROC_test,lateCycAVauROC_test);
        antiAnalysis.lateWinSI = cat(2,antiAnalysis.lateWinSI,lateWinSI);
        antiAnalysis.lateCycAV95CITest = cat(2,antiAnalysis.lateCycAV95CITest,lateCycSI95CITest');
        antiAnalysis.lateWinAV95CITest = cat(2,antiAnalysis.lateWinAV95CITest,lateWinSI95CITest');
        antiAnalysis.lateCycAVShuffTest = cat(2,antiAnalysis.lateCycAVShuffTest,lateCycAVShuffTest');
        antiAnalysis.lateWinAVShuffTest = cat(2,antiAnalysis.lateWinAVShuffTest,lateWinAVShuffTest');
        antiAnalysis.lateCycFracSI = cat(2,antiAnalysis.lateCycFracSI,lateCycFracSI);
        antiAnalysis.firstCycTC{visualTrials} = cat(2,antiAnalysis.firstCycTC{visualTrials},...
            mean(firstCycTC_vis,3));
        antiAnalysis.firstCycTC{auditoryTrials} = cat(2,antiAnalysis.firstCycTC{auditoryTrials},...
            mean(firstCycTC_aud,3));
        antiAnalysis.firstCycTC{allTrialsInd} = cat(2,antiAnalysis.firstCycTC{allTrialsInd},...
            mean(cat(3,firstCycTC_vis,firstCycTC_aud),3));
        antiAnalysis.firstCycTCErr{visualTrials} = cat(2,antiAnalysis.firstCycTCErr{visualTrials},...
            ste(firstCycTC_vis,3));
        antiAnalysis.firstCycTCErr{auditoryTrials} = cat(2,antiAnalysis.firstCycTCErr{auditoryTrials},...
            ste(firstCycTC_aud,3));
        antiAnalysis.firstCycTCErr{allTrialsInd} = cat(2,antiAnalysis.firstCycTCErr{allTrialsInd},...
            ste(cat(3,firstCycTC_vis,firstCycTC_aud),3));

        targetAnalysis.tc{1,visualTrials} = cat(2,targetAnalysis.tc{1,visualTrials},...
            mean(targetDataExpt(iexp).av(visualTrials).tc{1},3));
        targetAnalysis.tc{2,visualTrials} = cat(2,targetAnalysis.tc{2,visualTrials},...
            mean(targetDataExpt(iexp).av(visualTrials).tc{2},3));
        targetAnalysis.tc{3,visualTrials} = cat(2,targetAnalysis.tc{3,visualTrials},...
            mean(cat(3,targetDataExpt(iexp).av(visualTrials).tc{1},...
            targetDataExpt(iexp).av(visualTrials).tc{2}),3));
        targetAnalysis.tc{1,auditoryTrials} = cat(2,targetAnalysis.tc{1,auditoryTrials},...
            mean(targetDataExpt(iexp).av(auditoryTrials).tc{1},3));
        targetAnalysis.tc{2,auditoryTrials} = cat(2,targetAnalysis.tc{2,auditoryTrials},...
            mean(targetDataExpt(iexp).av(auditoryTrials).tc{2},3));
        targetAnalysis.tc{3,auditoryTrials} = cat(2,targetAnalysis.tc{3,auditoryTrials},...
            mean(cat(3,targetDataExpt(iexp).av(auditoryTrials).tc{1},...
            targetDataExpt(iexp).av(auditoryTrials).tc{2}),3));
        targetAnalysis.trOutTC(:,visualTrials) = cellfun(@(x,y) cat(2,x,mean(y(1:(nBaselineFr*2),:,:),3)),...
            targetAnalysis.trOutTC(:,visualTrials),...
            targetDataExpt(iexp).av(visualTrials).trOutTC','unif',0);
        targetAnalysis.trOutTC(:,auditoryTrials) = cellfun(@(x,y) cat(2,x,mean(y(1:(nBaselineFr*2),:,:),3)),...
            targetAnalysis.trOutTC(:,auditoryTrials),...
            targetDataExpt(iexp).av(auditoryTrials).trOutTC','unif',0);

        targetAnalysis.targets{visualTrials} = cat(2,targetAnalysis.targets{visualTrials},...
            cat(1,ones(1,nCells).*mean(targetDataExpt(iexp).av(visualTrials).targets{1}),...
            ones(1,nCells).*mean(targetDataExpt(iexp).av(visualTrials).targets{2})));
        targetAnalysis.targetsAuROC = cat(1,targetAnalysis.targetsAuROC,auroc_target);
        targetAnalysis.targetsAuROCTest = cat(1,targetAnalysis.targetsAuROCTest,auroc_target_test);
        taskTuningTest = cat(1,taskTuningTest,taskStimAnova' < 0.05);
        taskTuningPref = cat(1,taskTuningPref,taskTuningID);

        firstStimAuROC = cat(1,firstStimAuROC,auroc_first);
        lateStimAuROC = cat(1,lateStimAuROC,auroc_late);
        firstStimAuROCTest = cat(1,firstStimAuROCTest,auroc_first_test);
        lateStimAuROCTest = cat(1,lateStimAuROCTest,auroc_late_test);

        firstStimAuROC_aud = cat(1,firstStimAuROC_aud,auroc_first_aud);
        lateStimAuROC_aud = cat(1,lateStimAuROC_aud,auroc_late_aud);
        firstStimAuROCTest_aud = cat(1,firstStimAuROCTest_aud,auroc_first_aud_test);
        lateStimAuROCTest_aud = cat(1,lateStimAuROCTest_aud,auroc_late_aud_test);

    %     ind = taskStimAnova' > 0.05 & ttest(firstVisResp,tarVisResp{3},'dim',2);
    %     taskTuningID_temp = taskTuningID;
    %     taskTuningID_temp(ind) = 4;
    %     taskTuningPref_temp =cat(1,taskTuningPref_temp,taskTuningID_temp);
    end
    % antiAnalysis.adapt = cellfun(@(x,y) ...
    %     mean(x(respwin,:),1)./mean(y(respwin,:),1),...
    %     antiAnalysis.lateCycTC,antiAnalysis.firstCycTC,'unif',0);
    
    pctModCellsExpt = nModCellsPerExpt./cellfun(@sum,{respCellsExpt.lateCycRespCells});
    fprintf('%s+/-%s mod cells, rew mice\n',num2str(mean(pctModCellsExpt(logical(rewExptInd)))),...
        num2str(ste(pctModCellsExpt(logical(rewExptInd)),2)))
    fprintf('%s+/-%s mod cells, no rew mice\n',num2str(mean(pctModCellsExpt(logical(~rewExptInd)))),...
        num2str(ste(pctModCellsExpt(logical(~rewExptInd)),2)))

    cellInfo = struct;
    cellInfo.firstRespCells = logical(cell2mat({respCellsExpt.firstRespCells}'));
    cellInfo.lateRespCells = logical(cell2mat({respCellsExpt.lateRespCells}'));
    cellInfo.lateSuppCells = logical(cell2mat({respCellsExpt.lateSuppCells}'));
    cellInfo.lateCycRespCells = logical(cell2mat({respCellsExpt.lateCycRespCells}'));
    cellInfo.minRespCells = (mean(antiAnalysis.firstCycTC{visualTrials}(respwin,:),1) > ...
        minRespThreshold & mean(antiAnalysis.firstCycTC{auditoryTrials}(respwin,:),1) > ...
        minRespThreshold)';
    cellInfo.targetRespCells = logical(cell2mat({respCellsExpt.targetRespCells}'));
    cellInfo.isShortCycExpt = isShortCycExpt;
    cellInfo.isTuned = logical(cell2mat({oriTuningExpt.isTuned}))';
    cellInfo.oriResp = cell2mat({oriTuningExpt.oriResp}');
    cellInfo.oriRespErr = cell2mat({oriTuningExpt.oriRespErr}');
    cellInfo.oriFit = cell2mat({oriTuningExpt.fit})';
    cellInfo.oriPref = cell2mat({oriTuningExpt.oriPref})';
    cellInfo.hwhm = hwhmFromOriFit(cellInfo.oriFit(:,1:180)',1:180)';
    cellInfo.taskTuningPref = taskTuningPref;
    % cellInfo.taskTuningPref = taskTuningPref_temp;
    cellInfo.taskTuningTest = taskTuningTest;

    auroc_first = nan(length(cellInfo.firstRespCells),1);
    auroc_late = nan(length(cellInfo.firstRespCells),1);
    auroc_first_test = nan(length(cellInfo.firstRespCells),1);
    auroc_late_test = nan(length(cellInfo.firstRespCells),1);
    auroc_late_HT = lateStimAuROC(:,1);
    auroc_late_HT_test = lateStimAuROCTest(:,1) < 0.05;
    for icell = 1:length(cellInfo.firstRespCells)
        if ~cellInfo.targetRespCells(icell)
            auroc_first(icell) = firstStimAuROC(icell,3);
            auroc_late(icell) = lateStimAuROC(icell,3);
            auroc_first_test(icell) = firstStimAuROCTest(icell,3) < 0.05;
            auroc_late_test(icell) = lateStimAuROCTest(icell,3) < 0.05;
        elseif ~targetAnalysis.targetsAuROCTest(icell)
            auroc_first(icell) = firstStimAuROC(icell,3);
            auroc_late(icell) = lateStimAuROC(icell,3);
            auroc_first_test(icell) = firstStimAuROCTest(icell,3) < 0.05;
            auroc_late_test(icell) = lateStimAuROCTest(icell,3) < 0.05;
        elseif targetAnalysis.targetsAuROCTest(icell) && targetAnalysis.targetsAuROC(icell) < 0.5
            auroc_first(icell) = firstStimAuROC(icell,1);
            auroc_late(icell) = lateStimAuROC(icell,1);
            auroc_first_test(icell) = firstStimAuROCTest(icell,1) < 0.05;
            auroc_late_test(icell) = lateStimAuROCTest(icell,1) < 0.05;
        elseif targetAnalysis.targetsAuROCTest(icell) && targetAnalysis.targetsAuROC(icell) >= 0.5
            auroc_first(icell) = firstStimAuROC(icell,2);
            auroc_late(icell) = lateStimAuROC(icell,2);
            auroc_first_test(icell) = firstStimAuROCTest(icell,2) < 0.05;
            auroc_late_test(icell) = lateStimAuROCTest(icell,2) < 0.05;
        end
        
    end
    cellInfo.firstStimAuROC = auroc_first;
    cellInfo.lateStimAuROC = auroc_late;
    cellInfo.firstStimAuROCTest = auroc_first_test;
    cellInfo.lateStimAuROCTest = auroc_late_test;
    cellInfo.audFirstStimAuROC = firstStimAuROC_aud;
    cellInfo.audLateStimAuROC = lateStimAuROC_aud;
    cellInfo.audFirstStimAuROCTest = firstStimAuROCTest_aud < 0.05;
    cellInfo.audLateStimAuROCTest = lateStimAuROCTest_aud < 0.05;
    if strcmp(ds,'FSAV_attentionV1')
        cellInfo.dcModelCells = logical(cell2mat({decodeAnalysis.cellInd}'));
        cellInfo.targetWeight = targetWeightBoot;
        cellInfo.detectWeight = detectWeightBoot;
        cellInfo.targetComboWeight = nanmean(cell2mat(w_target_combo_boot'),2)';
        cellInfo.detectComboWeight = nanmean(cell2mat(w_detect_combo_boot'),2)';
        save([fnout 'imgAnalysisData'],...
            'decodeAnalysis','antiAnalysis','targetAnalysis','cellInfo')
    else
        cellInfo.targetWeight = targetWeightBoot;
        cellInfo.dcModelCells = logical(cell2mat({decodeAnalysis.cellInd}'));
    end
    save([fnout 'imgAnalysisData'],...
        'decodeAnalysis','antiAnalysis','targetAnalysis','cellInfo')
end
%% plotting params
respTCLim = [-0.005 0.05];
cycTCLim = [-0.005 0.015];
cycTCLim_minRespCells = [-0.005 0.025];
scatLim_win = [-0.2 0.6];
scatLim_cyc = [-0.035 0.085];
hmLim = [-0.1 0.1];
exCellTCLim = [-0.02 0.15];
oriRespLim = [-0.05 0.15];
siLim = [-10 10];
siOriLim = [-3 3];
oriBarLim_win = [0 0.08];
oriBarLim_resp = [0 0.04];
oriLim_taskResp = [-0.005 0.035];
oriNLim = [0 120];
oriTCLim = [-0.005 0.08];
targetTCLim = [-0.015 0.08];
outTCLim = [-0.005 0.04];
firstTCLim = [-0.005 0.04];
% adaptLim = [0 1];
suppTCLim = [-0.05 0.005];
suppScatLim_win = [-0.2 0.1];
suppScatLim_cyc = [-0.015 0.015];
cellRespTCLim = [-0.05 0.15];
exCellRespTCLim = [-0.1 0.1];
siBinLim = [-4 4];
stimRespLim = [-0.01 0.05];

tcStartFrame = 26;
cycTCEndTimeMs = 350;
cycTCEndFr = 45;
ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
ttLabel_target = -1000:250:900;
preTargetStimLabel = -700:350:0;
nFr_long = size(antiAnalysis.longTC{1,1},1);
tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_cyc = ((ttLabel_cyc./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr_target)+1);

nFr_cyc = size(antiAnalysis.lateCycTC{1,1},1);
tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
tt_targetTC = ((1:nFr_cyc)-(nBaselineFr+nVisDelayFr_target)).*(1000/frameRateHz);

lateWinTT = ([lateWinFr(1) lateWinFr(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT_target = (...
    [respwin_target(1) respwin_target(end)] - (nBaselineFr+nVisDelayFr_target))...
    .*(1000/frameRateHz);
baseWinTT = (...
    [basewin_0(1) basewin_0(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);

movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelMs = (movWinLabelFr - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);

weightLim = [-3 4.2];
binnedWeightLim = [-0.4 0.4];
weightLimSum = [-0.8 0.8];
siLimSum = [-0.5 2.5];


% siBins = [-10,-1,1,10];
siBins = [-10,0,10];
nSIBins = length(siBins)-1;
SI = antiAnalysis.lateCycSI;
[n,~,siBinID] = histcounts(SI,siBins);
% siBinID = nan(1,length(SI));
% siBinID(~antiAnalysis.lateCycAV95CITest) = 2;
% siBinID(antiAnalysis.lateCycAV95CITest & SI < 0) = 1;
% siBinID(antiAnalysis.lateCycAV95CITest & SI > 0) = 3;
% nSIBins = 3;

aurocBinID = nan(1,length(SI));
aurocBinID(~cellInfo.lateStimAuROCTest) = 3;
aurocBinID(cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC < 0.5) = 1;
aurocBinID(cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC > 0.5) = 2;


aurocBinID_aud = nan(1,length(SI));
aurocBinID_aud(~cellInfo.audLateStimAuROCTest) = 3;
aurocBinID_aud(cellInfo.audLateStimAuROCTest & cellInfo.audLateStimAuROC < 0.5) = 1;
aurocBinID_aud(cellInfo.audLateStimAuROCTest & cellInfo.audLateStimAuROC > 0.5) = 2;

Dr = cellfun(@(x) mean(x(respwin,:),1),antiAnalysis.lateCycTC,'unif',0);
Tr_vis = cellfun(@(x) mean(x(respwin_target,:),1),targetAnalysis.tc(:,visualTrials),'unif',0);
Tr_aud = cellfun(@(x) mean(x(respwin_target,:),1),targetAnalysis.tc(:,auditoryTrials),'unif',0);

%% plot anticipation analysis (Figure 2)

setFigParams4Print('landscape')

% heatmaps
lateWinRespAll = mean(antiAnalysis.longTC{allTrialsInd}(lateWinFr,:),1);
[~,lateWinSortInd] = sort(lateWinRespAll);
hm = flipud(antiAnalysis.longTC{allTrialsInd}(:,lateWinSortInd)');

figure
colormap(brewermap([],'*RdBu'));
% subplot 121
imagesc(hm(:,tcStartFrame:end))
hold on
if strcmp(ds,'FSAV_attentionV1')
    exCellInd = [exampleCell_1,exampleCell_2,exampleCell_3];
    exCellMat = zeros(1,length(cellInfo.firstRespCells));
    exCellMat(exCellInd) = 1;
    exCellSortInd = find(flip(exCellMat(lateWinSortInd)));
    hline(exCellSortInd,'k-')
end
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabelFr_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title('All Trials, All Cells')
print([fnout 'heatmapAllTrialsAllCells'],'-dpdf','-fillpage')

figure
colormap gray
% subplot 122
ind = cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells;
indBW = ind; indBW(ind) = 0; indBW(~ind) = 1;
imagesc(flipud(indBW(lateWinSortInd)))
figYAxis([],'Cell #',[])
figAxForm
colorbar
title('Anti. Resp. Cells')
print([fnout 'cellIDforHeatmap'],'-dpdf','-fillpage')

% time-courses and quantification scatters
setFigParams4Print('portrait')
figure
subplot 311
ind = cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells;
for iav = 1:2
    y = mean(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(ind) & ...
        cellInfo.isShortCycExpt),2);
    yerr = ste(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(ind) & cellInfo.isShortCycExpt),2);
    hold on
    shadedErrorBar_chooseColor(tt_longTC,y,yerr,cueColor{iav});
end
figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
figYAxis([],'dF/F',respTCLim)  
vline(lateWinTT,'k--')
hline(0,'k:')
figAxForm([],0)
title(sprintf('First Stim Responsive Cells (%s/%s)',...
    num2str(sum((...
    cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)...
    & cellInfo.isShortCycExpt)),...
    num2str(sum(cellInfo.isShortCycExpt))))

subplot 323
x = mean(antiAnalysis.longTC{visualTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
xerr = ste(x,2);
y = mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
yerr = ste(y,2);
plot(x,y,'.')
hold on
errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.')
plot(scatLim_win,scatLim_win,'k--')
plot(scatLim_win,[0 0],'k:')
plot([0 0],scatLim_win,'k:')
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',scatLim_win)
figYAxis([],'Auditory (dF/F)',scatLim_win)
figAxForm
title(sprintf('Late Window, All Resp. Cells (%s/%s), p = %s',...
    num2str(sum(cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),...
    num2str(length(cellInfo.firstRespCells)),...
    num2str(round(p,2,'significant'))))


subplot 325
for iav = 1:2
    y = mean(antiAnalysis.lateCycTC{iav}...
        ((tcStartFrame:end),(cellInfo.lateCycRespCells)),2);
    yerr = ste(antiAnalysis.lateCycTC{iav}...
        ((tcStartFrame:end),(cellInfo.lateCycRespCells)),2);
    hold on
    shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
end
figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
figYAxis([],'dF/F',cycTCLim)  
vline(lateWinTT,'k--')
hline(0,'k:')
vline(respWinTT,'k--')
figAxForm
title(sprintf('Late Cycle Resp. Cells (%s/%s)',...
    num2str(sum(cellInfo.lateCycRespCells)),...
    num2str(length(cellInfo.firstRespCells))))

subplot 326
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
xerr = ste(x,2);
y = mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
yerr = ste(y,2);
plot(x,y,'.')
hold on
errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.')
plot(scatLim_cyc,scatLim_cyc,'k--')
plot(scatLim_cyc,[0 0],'k:')
plot([0 0],scatLim_cyc,'k:')
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',scatLim_cyc)
figYAxis([],'Auditory (dF/F)',scatLim_cyc)
figAxForm
title(sprintf('Late Cycles, Late Cyc. Resp. Cells (%s/%s), p = %s',...
    num2str(sum(cellInfo.lateCycRespCells)),...
    num2str(length(cellInfo.firstRespCells)),...
    num2str(round(p,2,'significant'))))

print([fnout 'tcLongTrialsAndLateCycWithQuant_LateCycRespCells'],'-dpdf','-fillpage')

% example cells time-courses
if strcmp(ds,'FSAV_attentionV1')
    setFigParams4Print('landscape')
    exCellInd = [exampleCell_1,exampleCell_2,exampleCell_3];
    figure
    for icell = 1:3
        exampleCell = exCellInd(icell);
        if icell == 1
            plotOffset = 0;
        else
            plotOffset = plotOffset+4;
        end
        subplot(3,4,1+plotOffset)
        y = antiAnalysis.longTC{allTrialsInd}(tcStartFrame:end,exampleCell);
        yerr = antiAnalysis.longTCErr{allTrialsInd}(tcStartFrame:end,exampleCell);
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
        hold on
        figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
        if icell == 3
            figYAxis([],'dF/F',exCellTCLim-0.06)
        else
            figYAxis([],'dF/F',exCellTCLim)  
        end
        vline(lateWinTT,'k--')
        hline(0,'k:')
        figAxForm([],0)
        title(sprintf('Example Cell #%s',num2str(exampleCell)))
        
        subplot(3,4,2+plotOffset)
        y = antiAnalysis.firstCycTC{allTrialsInd}(tcStartFrame:end,exampleCell);
        yerr = antiAnalysis.firstCycTCErr{allTrialsInd}(tcStartFrame:end,exampleCell);
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,[0 0 0]);
        hold on
        figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
        if icell == 3
            figYAxis([],'dF/F',exCellTCLim-0.06)
        else
            figYAxis([],'dF/F',exCellTCLim)  
        end
        vline(lateWinTT,'k--')
        hline(0,'k:')
        figAxForm([],0)
        title(sprintf('Example Cell #%s',num2str(exampleCell)))
        
        subplot(3,4,3+plotOffset)
        y = antiAnalysis.lateCycTC{allTrialsInd}(tcStartFrame:end,exampleCell);
        yerr = antiAnalysis.lateCycTCErr{allTrialsInd}(tcStartFrame:end,exampleCell);
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,[0 0 0]);
        hold on
        figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
        if icell == 3
            figYAxis([],'dF/F',exCellTCLim-0.06)
        else
            figYAxis([],'dF/F',exCellTCLim)  
        end
        vline(lateWinTT,'k--')
        hline(0,'k:')
        figAxForm([],0)
        title(sprintf('Example Cell #%s',num2str(exampleCell)))

        subplot(3,4,4+plotOffset)
        y = cellInfo.oriResp(exampleCell,:);
        yerr = cellInfo.oriRespErr(exampleCell,:);
        errorbar(orientations,y,yerr,'.')
        hold on
        x = 0:180;
        y = cellInfo.oriFit(exampleCell,:);
        plot(x,y,'-')
        figXAxis([],'Orienation (deg)',[-10 190])
        figYAxis([],'dF/F',oriRespLim)
        figAxForm
        title(sprintf('Passive Ori. Tuning, Cell #%s',num2str(exampleCell)));
    end
    print([fnout 'exampleCellsTCWithTuning'],'-dpdf','-fillpage')
end

% selectivity
figure
subplot 121
y = antiAnalysis.lateCycSI(cellInfo.lateCycRespCells);
h = cdfplot(y);
hold on
vline(mean(y),'k-')
figXAxis([],'Selectivity Index',siLim)
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
title(sprintf('Late Cyc. Resp. Cells, mean = %s, ste = %s',...
    num2str(round(mean(y),2,'significant')),...
    num2str(round(ste(y,2),2,'significant'))))
subplot 122
y = antiAnalysis.lateCycSI_shuff(cellInfo.lateCycRespCells);
h = cdfplot(y);
hold on;
vline(mean(y),'k-')
figXAxis([],'Selectivity Index of Shuffle',siLim)
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
[~,p] = ttest(antiAnalysis.lateCycSI(cellInfo.lateCycRespCells), ...
    antiAnalysis.lateCycSI_shuff(cellInfo.lateCycRespCells));
title(sprintf('mean=%s, ste=%s, ttest p=%s',...
    num2str(round(mean(y),2,'significant')),...
    num2str(round(ste(y,2),2,'significant')),...
    num2str(round(p,2,'significant'))))

print([fnout 'selectivityCDF'],'-dpdf')

% VA auROC
figure
subplot 121
y = antiAnalysis.lateCycAVauROC(cellInfo.lateCycRespCells);
h = cdfplot(y);
hold on
vline(mean(y),'k-')
figXAxis([],'auROC_V_A',[0.3 0.7])
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
title(sprintf('Late Cyc. Resp. Cells, mean = %s, ste = %s',...
    num2str(round(mean(y),2,'significant')),...
    num2str(round(ste(y,2),2,'significant'))))

print([fnout 'VAauROCCDF'],'-dpdf')
%% create a structure with stats from analyses
imgStats = struct;
imgStats.nTrialsPerExpt = nTrialsPerExpt;
fprintf('%s +/- %s (%s-%s) trials per experiment\n',num2str(mean(nTrialsPerExpt)),...
    num2str(ste(nTrialsPerExpt,2)), num2str(min(nTrialsPerExpt)),...
    num2str(max(nTrialsPerExpt)))
imgStats.firstWinTimeMs = round(([firstWinFr(1) firstWinFr(end)]...
    - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz),4,'significant');
imgStats.lateWinTimeMs = round(lateWinTT,4,'significant');

%cell numbers
imgStats.nCells.totalCells = length(cellInfo.firstRespCells);
imgStats.nCells.allResp = sum(cellInfo.firstRespCells|cellInfo.lateCycRespCells|cellInfo.lateRespCells);
imgStats.nCells.firstResp = sum(cellInfo.firstRespCells);
imgStats.nCells.lateCycResp = sum(cellInfo.lateCycRespCells);
imgStats.nCells.lateOnlyResp = sum(...
    ~cellInfo.firstRespCells & (cellInfo.lateCycRespCells|cellInfo.lateRespCells));
% imgStats.nCells.lateTCResp = sum(cellInfo.lateRespCells);
imgStats.nCells.lateSupp = sum(cellInfo.lateSuppCells & ...
    ~cellInfo.lateCycRespCells & ~cellInfo.lateRespCells & ~cellInfo.firstRespCells);
imgStats.nCells.targetResp = sum(cellInfo.targetRespCells);
imgStats.nCells

x = mean(antiAnalysis.longTC{visualTrials}(firstWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
y = mean(antiAnalysis.longTC{auditoryTrials}(firstWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
imgStats.av(visualTrials).firstWinResp = mean(x);
imgStats.av(visualTrials).firstWinRespErr = ste(x,2);
imgStats.av(auditoryTrials).firstWinResp = mean(y);
imgStats.av(auditoryTrials).firstWinRespErr = ste(y,2);
[~,imgStats.firstWinTest] = ttest(x,y);
x = mean(antiAnalysis.longTC{visualTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
y = mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.lateCycRespCells)),1);
imgStats.av(visualTrials).lateWinResp = mean(x);
imgStats.av(visualTrials).lateWinRespErr = ste(x,2);
imgStats.av(auditoryTrials).lateWinResp = mean(y);
imgStats.av(auditoryTrials).lateWinRespErr = ste(y,2);
[~,imgStats.lateWinTest] = ttest(x,y);

% fprintf('Visual - Mean/Err First Win: %s/%s\n',...
%     num2str(round(imgStats.av(visualTrials).firstWinResp,2,'significant')),...
%     num2str(round(imgStats.av(visualTrials).firstWinRespErr,2,'significant')))
% fprintf('Auditory - Mean/Err First Win: %s/%s\n',...
%     num2str(round(imgStats.av(auditoryTrials).firstWinResp,2,'significant')),...
%     num2str(round(imgStats.av(auditoryTrials).firstWinRespErr,2,'significant')))
fprintf('First Win Test, p=%s\n',...
    num2str(round(imgStats.firstWinTest,2,'significant')))
% fprintf('Visual - Mean/Err Late Win: %s/%s\n',...
%     num2str(round(imgStats.av(visualTrials).lateWinResp,2,'significant')),...
%     num2str(round(imgStats.av(visualTrials).lateWinRespErr,2,'significant')))
% fprintf('Auditory - Mean/Err Late Win: %s/%s\n',...
%     num2str(round(imgStats.av(auditoryTrials).lateWinResp,2,'significant')),...
%     num2str(round(imgStats.av(auditoryTrials).lateWinRespErr,2,'significant')))
fprintf('Late Win Test, p=%s\n',...
    num2str(round(imgStats.lateWinTest,2,'significant')))

x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
y = mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,...
    (cellInfo.lateCycRespCells)),1);
imgStats.av(visualTrials).lateCycResp = mean(x);
imgStats.av(visualTrials).lateCycRespErr = ste(x,2);
imgStats.av(auditoryTrials).lateCycResp = mean(y);
imgStats.av(auditoryTrials).lateCycRespErr = ste(y,2);
[~,imgStats.lateCycTest] = ttest(x,y);

% fprintf('Visual - Mean/Err Late Cyc: %s/%s\n',...
%     num2str(round(imgStats.av(visualTrials).lateCycResp,2,'significant')),...
%     num2str(round(imgStats.av(visualTrials).lateCycRespErr,2,'significant')))
% fprintf('Auditory - Mean/Err Late Cyc: %s/%s\n',...
%     num2str(round(imgStats.av(auditoryTrials).lateCycResp,2,'significant')),...
%     num2str(round(imgStats.av(auditoryTrials).lateCycRespErr,2,'significant')))
fprintf('Late Cyc Test, p=%s\n',...
    num2str(round(imgStats.lateCycTest,2,'significant')))

imgStats.respCellsSI = antiAnalysis.lateCycSI(cellInfo.lateCycRespCells);
[~,imgStats.respCellsSITest] = ttest(antiAnalysis.lateCycSI(...
    cellInfo.lateCycRespCells));
% imgStats.nAVModCells = sum(antiAnalysis.lateCycAVTest(cellInfo.lateCycRespCells));
% fprintf('Mean/Err SI: %s/%s\n',...
%     num2str(round(mean(imgStats.respCellsSI,2),2,'significant')),...
%     num2str(round(ste(imgStats.respCellsSI,2),2,'significant')))
fprintf('SI Test, p=%s\n',...
    num2str(round(imgStats.respCellsSITest,2,'significant')))

if strcmp(ds,'FSAV_attentionV1')
    [~,imgStats.respCellsSITestVsNaive] = kstest2(imgStats.respCellsSI, ...
        imgStats_naive.respCellsSI,'tail','smaller');
    fprintf('Kolmogorov-Smirnov Test SI, behavior vs. naive, p = %s\n',...
        num2str(round(imgStats.respCellsSITestVsNaive,2,'significant')))
end
%% Figure 3 example cell responses    
iexp=10;
exCells = [13,17];
dmat = antiDataExpt(iexp).cycTC(visualTrials,5:end);
d_anti = [];
for i = 1:length(dmat)
    d_anti = cat(3,d_anti,dmat{i} - mean(dmat{i}(basewin_0,:,:),1));
end
d_tar = targetDataExpt(iexp).av(visualTrials).tc;
% ind = decodeAnalysis(iexp).cellInd;
% trOut = dc.outcome;
% stimID = discretize(dc.stim,oriBins);
stimLabel = {'D','HT','ET'};
figure
suptitle(sprintf('%s, example cell stim resp',antiDataExpt(iexp).exptName))
for icell = 1:2
    if icell == 1
        plotOffset = 0;
    else
        plotOffset = plotOffset+4;
    end 
    for istim = 1:3
        subplot(2,4,istim+plotOffset)
        if istim == 1
            r = squeeze(d_anti(:,exCells(icell),:));
            shadedErrorBar_chooseColor(tt_cycTC,mean(r(tcStartFrame:nFr_cyc,:),2),...
                ste(r(tcStartFrame:nFr_cyc,:),2),[0 0 0]);
            figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
        else
            r = squeeze(d_tar{istim-1}(:,exCells(icell),:));
            shadedErrorBar_chooseColor(tt_targetTC((tcStartFrame:nFr_cyc)-1),mean(r(tcStartFrame:nFr_cyc,:),2),...
                ste(r(tcStartFrame:nFr_cyc,:),2),[0 0 0]);
            figXAxis([],'Time (ms)',[tt_targetTC(tcStartFrame-1) cycTCEndTimeMs],[0,200],[0,200])
        end
%         trInd = stimID == istim;
%         r = squeeze(d_anti(:,exCells(icell),:));
        figYAxis([],'dF/F',exCellRespTCLim)  
        hline(0,'r:')
        vline(respWinTT,'k--')
        figAxForm
        title(sprintf('Cell %s, %s Resp.',num2str(exCells(icell)),...
            stimLabel{istim}))
    end
    subplot(2,4,4+plotOffset)
    x = orientations;
    y = oriTuningExpt(iexp).oriResp(exCells(icell),:);
    yerr = oriTuningExpt(iexp).oriRespErr(exCells(icell),:);
    errorbar(x,y,yerr,'.')
    hold on
    y = oriTuningExpt(iexp).fit(:,exCells(icell));
    plot(0:180,y,'-')
    title(sprintf('pref=%s, peak=%s, theta90=%s',...
        num2str(orientations(oriTuningExpt(iexp).oriPref(exCells(icell)))),...
        num2str(oriTuningExpt(iexp).fitPeak(exCells(icell))),...
        num2str(oriTuningExpt(iexp).tuningReliability(exCells(icell)))))
    figXAxis([],'Orientation',[-10 190],x,x)
    figYAxis([],'dF/F',[])
    figAxForm
end
print([fnout 'exCellResponses' decodeDataExpt(iexp).exptName],...
    '-dpdf','-fillpage')

%% plot orientation analysis (Figure 3)
if strcmp(ds,'FSAV_attentionV1')
    minCellN_SIFRmatch = 41;
else
    minCellN_SIFRmatch = 33;
end
% adaptAnalysis = struct;
% ind = cellInfo.firstRespCells & cellInfo.minRespCells;
% adaptAnalysis.firstRespCells = cellfun(@(x) mean(x(ind)),antiAnalysis.adapt);
% adaptAnalysis.firstRespCellsErr = cellfun(@(x) ste(x(ind),2),antiAnalysis.adapt);
% [~,adaptAnalysis.testAV] = ttest(antiAnalysis.adapt{visualTrials}(ind),...
%     antiAnalysis.adapt{auditoryTrials}(ind),'tail','right');
% ind2 = ind & cellInfo.isTuned;
% adaptAnalysis.oriGroupsTest = anova1(antiAnalysis.adapt{3}(ind2),...
%     cellInfo.oriPref(ind2),'off');

% get first, late win, and late cycle response for passive ori tuning groups
oriGroups = struct;
oriGroups.n = nan(1,nOri);
oriGroups.nShortCyc = nan(1,nOri);
oriGroups.nTarOrDist = nan(1,nOri);
oriGroups.nTar = nan(1,nOri);
oriGroups.lateCycSI = nan(1,nOri);
oriGroups.lateCycSIErr = nan(1,nOri);
oriGroups.lateCycAVauROC = nan(1,nOri);
oriGroups.lateCycAVauROCErr = nan(1,nOri);
oriGroups.lateCycAVauROCData = cell(1,nOri);
oriGroups.lateCycRespAll = nan(1,nOri);
oriGroups.lateCycRespAllErr = nan(1,nOri);
oriGroups.lateCycRespAllData = cell(1,nOri);
oriGroups.firstResp = nan(2,nOri);
oriGroups.firstRespErr = nan(2,nOri);
oriGroups.lateWin = nan(2,nOri);
oriGroups.lateWinErr = nan(2,nOri);
oriGroups.lateCycResp = nan(2,nOri);
oriGroups.lateCycRespErr = nan(2,nOri);
oriGroups.lateCycRespDiff = nan(1,nOri);
oriGroups.lateCycRespDiffErr = nan(1,nOri);
oriGroups.firstTC = cell(2,nOri);
oriGroups.firstTCErr = cell(2,nOri);
oriGroups.longTC = cell(2,nOri);
oriGroups.longTCErr = cell(2,nOri);
oriGroups.cycTC = cell(2,nOri);
oriGroups.cycTCErr = cell(2,nOri);
oriGroups.lateCycSITestEaOri = nan(1,nOri);
oriGroups.lateCycTest = nan(1,nOri);
oriGroups.lateWinTest = nan(1,nOri);
oriGroups.targetTC = cell(1,nOri);
oriGroups.targetTCErr = cell(1,nOri);
oriGroups.targetTuningResp = cell(1,nOri);
oriGroups.targetTuningRespErr = cell(1,nOri);
oriGroups.targetTuningStim = cell(1,nOri);
oriGroups.targetTuningStimErr = cell(1,nOri);
oriGroups.firstStimRespForTargetAnalysis = nan(1,nOri);
oriGroups.firstStimRespErrForTargetAnalysis = nan(1,nOri);
oriGroups.matchSIxFR = nan(1,nOri);
oriGroups.matchSIxFRErr = nan(1,nOri);
oriGroups.matchedFR = nan(1,nOri);
oriGroups.matchedFRErr = nan(1,nOri);
oriGroups.matchSIxFRTestData = cell(1,nOri);
oriGroups.taskTuningData = cell(nOri,3);
% adaptAnalysis.oriGroupsN = nan(1,nOri);
% adaptAnalysis.oriGroups = nan(3,nOri);
% adaptAnalysis.oriGroupsErr = nan(3,nOri);
for iori = 1:nOri
    ind = cellInfo.isTuned & cellInfo.oriPref == iori &...
        (cellInfo.lateCycRespCells);
%     ind2 = cellInfo.isTuned & cellInfo.oriPref == iori & ...
%         cellInfo.firstRespCells & cellInfo.minRespCells;
    ind_tarAndDist = cellInfo.isTuned & cellInfo.oriPref == iori &...
        (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    ind_tar = cellInfo.isTuned & cellInfo.oriPref == iori & cellInfo.targetRespCells;
    oriGroups.n(iori) = sum(ind);
    oriGroups.nShortCyc(iori) = sum(ind & cellInfo.isShortCycExpt);
    oriGroups.nTarOrDist(iori) = sum(ind_tarAndDist);
    oriGroups.nTar(iori) = sum(ind_tar);
%     adaptAnalysis.oriGroupsN(iori) = sum(ind2);
    oriGroups.lateCycSI(iori) = mean(antiAnalysis.lateCycSI(ind));
    oriGroups.lateCycSIErr(iori) = ste(antiAnalysis.lateCycSI(ind),2);
    oriGroups.lateCycAVauROC(iori) = mean(antiAnalysis.lateCycAVauROC(ind),2);
    oriGroups.lateCycAVauROCErr(iori) = ste(antiAnalysis.lateCycAVauROC(ind),2);
    oriGroups.lateCycAVauROCData{iori} = antiAnalysis.lateCycAVauROC(ind);
    [~,oriGroups.lateCycSITestEaOri(iori)] = ttest(antiAnalysis.lateCycSI(ind));
    [~,oriGroups.lateWinTest(iori)] = ttest(mean(antiAnalysis.longTC{visualTrials}...
        (lateWinFr,ind & cellInfo.isShortCycExpt),1),...
        mean(antiAnalysis.longTC{auditoryTrials}...
        (lateWinFr,ind & cellInfo.isShortCycExpt),1));
    [~,oriGroups.lateCycTest(iori)] = ttest(mean(antiAnalysis.lateCycTC{visualTrials}...
        (respwin,ind),1),mean(antiAnalysis.lateCycTC{auditoryTrials}...
        (respwin,ind),1));
    oriGroups.lateCycRespDiff(iori) = mean(mean(antiAnalysis.lateCycTC{visualTrials}...
        (respwin,ind),1) - mean(antiAnalysis.lateCycTC{auditoryTrials}...
        (respwin,ind),1),2);
    oriGroups.lateCycRespDiffErr(iori) = ste(mean(antiAnalysis.lateCycTC{visualTrials}...
        (respwin,ind),1) - mean(antiAnalysis.lateCycTC{auditoryTrials}...
        (respwin,ind),1),2);
    
    oriGroups.lateCycRespAll(iori) = mean(mean(...
            antiAnalysis.lateCycTC{3}(respwin,ind),1));
    oriGroups.lateCycRespAllErr(iori) = ste(mean(...
            antiAnalysis.lateCycTC{3}(respwin,ind),1),2);
    oriGroups.lateCycRespAllData{iori} = mean(...
            antiAnalysis.lateCycTC{3}(respwin,ind),1);
    
    oriGroups.targetTC{iori} = cell2mat(cellfun(@(x) mean(x(:,ind_tar),2),...
        targetAnalysis.tc(1:2,visualTrials),'unif',0)');
    oriGroups.targetTCErr{iori} = cell2mat(cellfun(@(x) ste(x(:,ind_tar),2),...
        targetAnalysis.tc(1:2,visualTrials),'unif',0)');
    
%     oriGroups.targetTuningResp{iori} = cat(1,mean(mean(...
%             antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind_tarAndDist),1)),cellfun(@(x) ...
%         mean(mean(x(respwin_target,ind_tarAndDist),1) - mean(x(basewin_0,ind_tarAndDist),1)),...
%         targetAnalysis.tc(1:2,visualTrials)));
%     oriGroups.targetTuningRespErr{iori} = cat(1,ste(mean(...
%             antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind_tarAndDist),1),2),cellfun(@(x) ...
%         ste(mean(x(respwin_target,ind_tarAndDist),1) - mean(x(basewin_0,ind_tarAndDist),1),2),...
%         targetAnalysis.tc(1:2,visualTrials)));
%     oriGroups.firstStimRespForTargetAnalysis(iori) = mean(mean(...
%             antiAnalysis.firstCycTC{visualTrials}(respwin,ind_tarAndDist),1));
%     oriGroups.firstStimRespErrForTargetAnalysis(iori) = ste(mean(...
%             antiAnalysis.firstCycTC{visualTrials}(respwin,ind_tarAndDist),1),2);
%     oriGroups.targetTuningStim{iori} = cat(1,0,...
%         mean(targetAnalysis.targets{visualTrials}(:,ind_tarAndDist),2));
%     oriGroups.targetTuningStimErr{iori} = cat(1,0,...
%         ste(targetAnalysis.targets{visualTrials}(:,ind_tarAndDist),2));
    
    oriGroups.taskTuningData{iori,1} = mean(...
            antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind),1)';
    oriGroups.taskTuningData(iori,2:3) = cellfun(@(x) ...
        mean(x(respwin_target,ind),1)' - mean(x(basewin_0,ind),1)',...
        targetAnalysis.tc(1:2,visualTrials),'unif',0);
    oriGroups.targetTuningResp{iori} = cat(1,mean(mean(...
            antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind),1)),cellfun(@(x) ...
        mean(mean(x(respwin_target,ind),1) - mean(x(basewin_0,ind),1)),...
        targetAnalysis.tc(1:2,visualTrials)));
    oriGroups.targetTuningRespErr{iori} = cat(1,ste(mean(...
            antiAnalysis.lateCycTC{visualTrials}(respwin_target,ind),1),2),cellfun(@(x) ...
        ste(mean(x(respwin_target,ind),1) - mean(x(basewin_0,ind),1),2),...
        targetAnalysis.tc(1:2,visualTrials)));
    oriGroups.firstStimRespForTargetAnalysis(iori) = mean(mean(...
            antiAnalysis.firstCycTC{visualTrials}(respwin,ind),1));
    oriGroups.firstStimRespErrForTargetAnalysis(iori) = ste(mean(...
            antiAnalysis.firstCycTC{visualTrials}(respwin,ind),1),2);
    oriGroups.targetTuningStim{iori} = cat(1,0,...
        mean(targetAnalysis.targets{visualTrials}(:,ind),2));
    oriGroups.targetTuningStimErr{iori} = cat(1,0,...
        ste(targetAnalysis.targets{visualTrials}(:,ind),2));
    
    lateCycResp_grp = Dr{3}(ind);
    lateCycSI_grp = antiAnalysis.lateCycSI(ind);
    [sortFR,sortFRInd] = sort(lateCycResp_grp,'descend');
    sortSI = lateCycSI_grp(sortFRInd);
    oriGroups.matchSIxFR(iori) = mean(sortSI(1:minCellN_SIFRmatch),2);
    oriGroups.matchSIxFRErr(iori) = ste(sortSI(1:minCellN_SIFRmatch),2);
    oriGroups.matchedFR(iori) = mean(sortFR(1:minCellN_SIFRmatch),2);
    oriGroups.matchedFRErr(iori) = ste(sortFR(1:minCellN_SIFRmatch),2);
    oriGroups.matchSIxFRTestData(iori) = {sortSI(1:minCellN_SIFRmatch)};
            
    for iav = 1:2
        oriGroups.firstResp(iav,iori) = mean(mean(...
            antiAnalysis.firstCycTC{iav}(respwin,ind),1));
        oriGroups.firstRespErr(iav,iori) = ste(mean(...
            antiAnalysis.firstCycTC{iav}(respwin,ind),1),2);
        
        oriGroups.firstTC(iav,iori) = {mean(...
            antiAnalysis.firstCycTC{iav}(:,ind),2)};
        oriGroups.firstTCErr(iav,iori) = {ste(...
            antiAnalysis.firstCycTC{iav}(:,ind),2)};
        
        oriGroups.lateWin(iav,iori) = mean(mean(...
            antiAnalysis.longTC{iav}(lateWinFr,ind),1));
        oriGroups.lateWinErr(iav,iori) = ste(mean(...
            antiAnalysis.longTC{iav}(lateWinFr,ind),1),2);
        
        oriGroups.lateCycResp(iav,iori) = mean(mean(...
            antiAnalysis.lateCycTC{iav}(respwin,ind),1));
        oriGroups.lateCycRespErr(iav,iori) = ste(mean(...
            antiAnalysis.lateCycTC{iav}(respwin,ind),1),2);
        
        oriGroups.longTC{iav,iori} = mean(...
            antiAnalysis.longTC{iav}(:,ind & cellInfo.isShortCycExpt),2);
        oriGroups.longTCErr{iav,iori} = ste(...
            antiAnalysis.longTC{iav}(:,ind & cellInfo.isShortCycExpt),2);
        
        oriGroups.cycTC{iav,iori} = mean(...
            antiAnalysis.lateCycTC{iav}(:,ind),2);
        oriGroups.cycTCErr{iav,iori} = ste(...
            antiAnalysis.lateCycTC{iav}(:,ind),2);
%         adaptAnalysis.oriGroups(iav,iori) = mean(...
%             antiAnalysis.adapt{iav}(ind2),2);
%         adaptAnalysis.oriGroupsErr(iav,iori) = ste(...
%             antiAnalysis.adapt{iav}(ind2),2);
%         if iav == 2
%             adaptAnalysis.oriGroups(3,iori) = mean(...
%                 antiAnalysis.adapt{3}(ind2),2);
%             adaptAnalysis.oriGroupsErr(3,iori) = ste(...
%                 antiAnalysis.adapt{3}(ind2),2);
%         end
        
        
    end
end
ind = cellInfo.isTuned &...
    (cellInfo.lateCycRespCells);
lateCycRespDiffData = ...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1) - ...
    mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,ind),1);
[oriGroups.lateCycRespTest,~,stats] = anova1(lateCycRespDiffData,cellInfo.oriPref(ind),'off');
oriGroups.lateCycRespTestPostHoc = multcompare(stats,[],'off');
% oriGroups.lateCycSITest = anova1(antiAnalysis.lateCycSI(ind),cellInfo.oriPref(ind),'off');
[lateCycSITest,~,stats] = anova1(antiAnalysis.lateCycSI(ind),cellInfo.oriPref(ind),'off');
% lmTest = fitlm(antiAnalysis.lateCycSI(ind),cellInfo.oriPref(ind));
oriGroups.lateCycSITest = lateCycSITest;
[~,oriGroups.FRmatchedSITestEaOri] = cellfun(@ttest,oriGroups.matchSIxFRTestData);

[lateCycAVauROCTest,~,stats] = anova1(antiAnalysis.lateCycAVauROC(ind),cellInfo.oriPref(ind),'off');
% lmTest = fitlm(antiAnalysis.lateCycSI(ind),cellInfo.oriPref(ind));
oriGroups.lateCycAVauROCTest = lateCycAVauROCTest;
[~,oriGroups.lateCycAVauROCTestEaOri] = cellfun(@(x) ttest(x,0.5),oriGroups.lateCycAVauROCData);

lateCycRespAV = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1);
[oriGroups.lateCycRespAllTest,~,stats] = anova1(lateCycRespAV,cellInfo.oriPref(ind),'off');
oriGroups.lateCycRespTestPostHoc = multcompare(stats);

ind = cellInfo.isTuned & (cellInfo.firstRespCells|cellInfo.lateCycRespCells | cellInfo.targetRespCells);
popTuningData = reshape(cat(1,mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),...
    cell2mat(cellfun(@(x) mean(x(respwin_target,ind),1) - mean(x(basewin_0,ind),1),...
        targetAnalysis.tc(1:2,visualTrials),'unif',0))),[3*sum(ind),1]);
popTuningGrp1 = repmat([1:3]',[sum(ind),1]);
popTuningGrp2 = reshape(repmat(cellInfo.oriPref(ind)',[3,1]),[3*sum(ind),1]);
oriGroups.popTuningTest = anovan(popTuningData,{popTuningGrp1,popTuningGrp2});

taskTuningTest = twoWayAnovaUnmatched(oriGroups.taskTuningData,{'Ori';'Stim'});
[p,~,stats] = anova1(Tr_vis{1}(cellInfo.lateCycRespCells),...
    cellInfo.oriPref(cellInfo.lateCycRespCells),'off');
tbl = multcompare(stats,[],'off');

% grpID = [];
% for igrp = 1:4
%     grpID = cat(2,grpID,ones(1,length(oriGroups.matchSIxFRTestData{igrp})).*igrp);
% end
% [oriGroups.matchSIxFRTest,~,stats] = anova1(cell2mat(...
%     oriGroups.matchSIxFRTestData),grpID,'off');

figure
suptitle('Tuned, First or Late Distractor Responsive Neurons')
for iori = 1:nOri
    subplot(2,2,iori)
    for iav = 1:2
        y = oriGroups.cycTC{iav,iori}(tcStartFrame:end);
        yerr = oriGroups.cycTCErr{iav,iori}(tcStartFrame:end);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',cycTCLim)  
    vline(respWinTT,'k--')
    hline(0,'k:')
    figAxForm
    title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
        num2str(oriGroups.n(iori)),num2str(round(oriGroups.lateCycTest(iori),2,'significant'))))
    
    
end

print([fnout 'tuningAnalysis_oriGroupsLateCycTC'],'-dpdf','-fillpage')    


figure
suptitle('First Cyc; Tuned, First or Late Distractor Responsive Neurons')
for iori = 1:nOri
    subplot(2,2,iori)
    for iav = 1:2
        y = oriGroups.firstTC{iav,iori}(tcStartFrame:end);
        yerr = oriGroups.firstTCErr{iav,iori}(tcStartFrame:end);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',firstTCLim)  
    vline(respWinTT,'k--')
    hline(0,'k:')
    figAxForm
    title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
        num2str(oriGroups.n(iori)),num2str(round(oriGroups.lateCycTest(iori),2,'significant'))))
    
end
print([fnout 'tuningAnalysis_oriGroupsFirstTC'],'-dpdf','-fillpage')  



setFigParams4Print('portrait')
figure;
suptitle('Tuned, First or Late Distractor Responsive Neurons, Short Cyc Expts')
for iori = 1:nOri
    subplot(4,2,iori)
    for iav = 1:2
        y = oriGroups.longTC{iav,iori}(tcStartFrame:end);
        yerr = oriGroups.longTCErr{iav,iori}(tcStartFrame:end);
        hold on
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,cueColor{iav});
    end
    figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
    figYAxis([],'dF/F',oriTCLim)  
    vline(lateWinTT,'k--')
    hline(0,'k:')
    figAxForm([],0)
    title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
        num2str(oriGroups.nShortCyc(iori)),num2str(round(oriGroups.lateWinTest(iori),2,'significant'))))
    
%     subplot(4,2,iori)
%     for iav = 1:2
%         y = oriGroups.cycTC{iav,iori}(tcStartFrame:end);
%         yerr = oriGroups.cycTCErr{iav,iori}(tcStartFrame:end);
%         hold on
%         shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
%     end
%     figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
%     figYAxis([],'dF/F',oriTCLim)  
%     vline(respWinTT,'k--')
%     hline(0,'k:')
%     figAxForm
%     title(sprintf('Pref %s, n=%s, p=%s',num2str(orientations(iori)),...
%         num2str(oriGroups.n(iori)),num2str(round(oriGroups.lateCycTest(iori),2,'significant'))))
    
    subplot(4,2,iori+4)
    for itar = 1:2
        y = oriGroups.targetTC{iori}(tcStartFrame:end,itar);
        yerr = oriGroups.targetTCErr{iori}(tcStartFrame:end,itar);
        hold on
        shadedErrorBar_chooseColor(tt_cycTC,y,yerr,hiLoColor{itar});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',targetTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    figAxForm
    title(sprintf('Pref %s, n=%s',num2str(orientations(iori)),...
        num2str(oriGroups.nTar(iori))))
    
end

print([fnout 'tuningAnalysis_oriGroupsTC'],'-dpdf','-fillpage')    

setFigParams4Print('portrait')
figure;
suptitle(sprintf('Tuned, Late Distractor-Responsive Neurons, %s/%s',...
    num2str(sum(cellInfo.isTuned & (cellInfo.lateCycRespCells))),...
    num2str(length(cellInfo.isTuned))))
xsub = [-0.25 +0.25];
subplot 231
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.firstResp(iav,:),...
        oriGroups.firstRespErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title('First Stim')

subplot 232
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.lateWin(iav,:),...
        oriGroups.lateWinErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_win)
figAxForm
title('Late Window')

subplot 233
for iav = 1:2
    hold on
    errorbar((1:nOri)+xsub(iav),oriGroups.lateCycResp(iav,:),...
        oriGroups.lateCycRespErr(iav,:),'.')
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriBarLim_resp)
figAxForm
title(sprintf('Late Stim Resp, One-Way ANOVA p=%s',num2str(round(...
    oriGroups.lateCycRespTest,2,'significant'))))

subplot 234
h = bar(oriGroups.n,'EdgeColor','none','BarWidth',0.5);
hold on
for iori = 1:nOri
    text(iori,oriGroups.n(iori)+1,num2str(oriGroups.n(iori)))
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'dF/F',oriNLim)
figAxForm

subplot 235
h = bar(oriGroups.lateCycSI,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:nOri,oriGroups.lateCycSI,oriGroups.lateCycSIErr,'.')
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'Selectivity',siOriLim)
figAxForm
title(sprintf('Late Stim Resp, One-Way ANOVA p=%s',num2str(round(...
    oriGroups.lateCycSITest,2,'significant'))))

subplot 236
% h = bar(oriGroups.lateCycAVauROC,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:nOri,oriGroups.lateCycAVauROC,oriGroups.lateCycAVauROCErr,'.')
for iori = 1:4
    text(iori,0.55,num2str(oriGroups.lateCycAVauROCTestEaOri(iori)))
end
figXAxis([],'Pref. Ori. (deg)',[0 nOri+1],[1:nOri],orientations)
figYAxis([],'auROC_V_A',[0.45 0.55])
figAxForm
hline(0.5,'r:')
title(sprintf('Late Stim Resp, One-Way ANOVA p=%s',num2str(round(...
    oriGroups.lateCycAVauROCTest,2,'significant'))))

print([fnout 'tuningAnalysis_lateDistRespCells'],'-dpdf','-fillpage')

% population task tuning curves of ori pref groups
setFigParams4Print('landscape')    
figure
suptitle({sprintf('Tuned and Distractor or Target Resp. Neurons %s/%s',...
    num2str(sum(cellInfo.isTuned & ...
    (cellInfo.firstRespCells | cellInfo.lateRespCells | cellInfo.targetRespCells))),...
    num2str(length(cellInfo.isTuned)));sprintf(...
    '2-Way ANOVA: Main Effect Task Stim, p=%s; Main Effect Ori Pref, p=%s',...
    num2str(oriGroups.popTuningTest(1)),num2str(oriGroups.popTuningTest(2)))})
for iori = 1:nOri
    subplot(2,2,iori)
    x = oriGroups.targetTuningStim{iori};
    xerr = oriGroups.targetTuningStimErr{iori};
    y = oriGroups.targetTuningResp{iori};
    yerr = oriGroups.targetTuningRespErr{iori};
    errorbar(0,oriGroups.firstStimRespForTargetAnalysis(iori),...
        oriGroups.firstStimRespErrForTargetAnalysis(iori),'.')
    hold on
    errorbar(x,y,yerr,yerr,xerr,xerr,'.')
    figXAxis([],'Task Ori. (deg)',[-10 100],x,oriBins([1,3,4]))
    figYAxis([],'dF/F',oriLim_taskResp)  
    vline(lateWinTT,'k--')
    hline(0,'k:')
    figAxForm([],0)
    title(sprintf('Pref %s, n=%s',num2str(orientations(iori)),...
        num2str(oriGroups.n(iori))))
end

print([fnout 'tuningAnalysis_oriGroupsTaskTuning'],'-dpdf','-fillpage')   

% orientation groups binned by SI
ind = cellInfo.lateCycRespCells & cellInfo.isTuned;
frOriInBin = nan(nSIBins,4);
frSIInOriBin = nan(4,nSIBins);
siInOriBinChiTest = nan(4,1);
% binnedSIxOri = nan(4,2);
for ibin = 1:nSIBins
    binInd = siBinID' == ibin & ind;
    frOriInBin(ibin,:) = histcounts(cellInfo.oriPref(binInd),4,...
        'Normalization','probability');
%     frOriInBin(ibin,:) = histcounts(cellInfo.oriPref(binInd),4);
end
for iori = 1:4
    binInd = ind & cellInfo.oriPref == iori;
    frSIInOriBin(iori,:) = histcounts(siBinID(binInd),1:(nSIBins+1),...
        'Normalization','probability');
%     binnedSIxOri{iori,1} = SI(binInd & SI' < 0);
    n = sum(binInd);
    n2 = sum(binInd' & siBinID == 1);
    [~,siInOriBinChiTest(iori)] = prop_test([n2, round(n/2)],[n, n], 0.05);
end
binnedHT = cell(4,nSIBins);
binnedET = cell(4,nSIBins);
binnedT = cell(4,nSIBins);
binnedD = cell(4,nSIBins);
binnedSI = cell(4,nSIBins);
for ibin = 1:nSIBins
    for iori = 1:4
        binInd = ind & siBinID' == ibin & cellInfo.oriPref == iori;
        binnedHT{iori,ibin} = Tr_vis{1}(binInd);
        binnedET{iori,ibin} = Tr_vis{2}(binInd);
        binnedT{iori,ibin} = Tr_vis{3}(binInd);
        binnedD{iori,ibin} = Dr{3}(binInd);
        binnedSI{iori,ibin} = SI(binInd);
    end
end

% matName = {'Stim';'SI'};
oriSIBinTest = nan(1,nSIBins);
for ibin = 1:nSIBins
    binInd = ind & siBinID' == ibin;
    oriSIBinTest(ibin)= anova1(SI(binInd),cellInfo.oriPref(binInd),'off');
end
% for iori = 1:4
%     d = cat(1,binnedD(iori,:),binnedT(iori,:));
%     oriSIBinTest{iori} = twoWayAnovaUnmatched(d,matName);
% end

for iori = 1:4
    fprintf('%s deg - SI x D/T Two-Way ANOVA, p: %s=%s, %s=%s int=%s\n',...
        num2str(orientations(iori)),...
        oriSIBinTest{iori}.varTestName{1},num2str(round(oriSIBinTest{iori}.p(1),2,'significant')),...
        oriSIBinTest{iori}.varTestName{2},num2str(round(oriSIBinTest{iori}.p(2),2,'significant')),...
        num2str(round(oriSIBinTest{iori}.p(3),2,'significant')))
    if any(oriSIBinTest{iori}.p < 0.05)
%         disp(oriSIBinTest{iori}.posthoc)
        for ip = 1:2
            disp(oriSIBinTest{iori}.posthoc{ip})
        end
    end
end



figure
subplot 221
bar(frOriInBin,'stacked')
figXAxis([],'',[0 4],1:2,{'-SI','+SI'})
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
legend(cellfun(@num2str,num2cell(orientations),'unif',0),...
    'location','northeastoutside')
subplot 222
bar(frSIInOriBin,'stacked')
oriLabel = cellfun(@num2str,num2cell(orientations),'unif',0);
figXAxis([],'',[0 5],1:4,oriLabel)
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
legend({'-SI','+SI'},...
    'location','northeastoutside')
subplot 223
for ibin = 1:nSIBins
    y = cellfun(@mean,binnedSI(:,ibin));
    yerr = cellfun(@(x) ste(x,2),binnedSI(:,ibin));
    hold on
    errorbar(1:4,y,yerr,'.')
end
figXAxis([],'',[0 5],1:4,oriLabel)
figYAxis([],'SI',siBinLim)
figAxForm
legend({'-SI','+SI'},...
    'location','northeastoutside')

% for iori = 1:4
%     subplot(3,2,iori+2)
%     x = cellfun(@mean,binnedSI(iori,:));
%     xerr = cellfun(@(x) ste(x,2),binnedSI(iori,:));
%     y = cellfun(@mean,binnedD(iori,:));
%     yerr = cellfun(@(x) ste(x,2),binnedD(iori,:));    
%     hold on
%     errorbar(x,y,yerr,yerr,xerr,xerr,'.-')      
% %     y = cellfun(@mean,binnedHT(iori,:));
% %     yerr = cellfun(@(x) ste(x,2),binnedHT(iori,:));    
% %     errorbar(x,y,yerr,yerr,xerr,xerr,'.-')      
% %     y = cellfun(@mean,binnedET(iori,:));
% %     yerr = cellfun(@(x) ste(x,2),binnedET(iori,:));    
% %     errorbar(x,y,yerr,yerr,xerr,xerr,'.-')    
%     y = cellfun(@mean,binnedT(iori,:));
%     yerr = cellfun(@(x) ste(x,2),binnedT(iori,:));    
%     errorbar(x,y,yerr,yerr,xerr,xerr,'.-')  
%     title(sprintf('Pref. %s', num2str(orientations(iori))))
%     figXAxis([],'SI (mean in bin)',siBinLim)
%     figYAxis([],'Stim Resp (dF/F)',stimRespLim)
%     figAxForm
%     L=legend({'D_a_l_l','T_a_l_l'},'location','northeastoutside');
% end
print([fnout 'tuningAnalysis_oriGroupsSIBinning'],'-dpdf','-fillpage')  

% %SI matched for FR
% figure
% suptitle(sprintf('Use top %s cells from each group (ranked by dF/F resp to late vis stim',...
%     num2str(minCellN_SIFRmatch)))
% subplot 121
% errorbar(oriGroups.matchSIxFR,oriGroups.matchSIxFRErr,'.')
% hold on
% figXAxis([],'Ori. Pref.',[0 nOri+1],1:nOri,orientations)
% figYAxis([],'Matched Selectivity',[-1 3])
% hline(0,'k:')
% figAxForm
% subplot 122
% errorbar(oriGroups.matchedFR,oriGroups.matchedFRErr,'.')
% hold on
% figXAxis([],'Ori. Pref.',[0 nOri+1],1:nOri,orientations)
% figYAxis([],'Matched dF/F',[])
% hline(0,'k:')
% figAxForm
% print([fnout 'tuningAnalysis_oriGroupsMatchedSIxFR'],'-dpdf','-fillpage')   

%%
imgStats.nCells.oriGroups = oriGroups.n;
fprintf('n 0/45/90/135 pref = %s\n', num2str(imgStats.nCells.oriGroups))

imgStats.lateCycRespAllxOri = oriGroups.lateCycRespAll;
imgStats.lateCycRespAllxOriErr = oriGroups.lateCycRespAllErr;
% fprintf('Late Cyc Resp (All Trials) 0/45/90/135 pref = %s\n', ...
%     num2str(round(imgStats.lateCycRespAllxOri,2,'significant')))
% fprintf('Late Cyc Resp Err (All Trials) 0/45/90/135 pref = %s\n', ...
%     num2str(round(imgStats.lateCycRespAllxOriErr,2,'significant')))

imgStats.lateCycRespDiffxOri = oriGroups.lateCycRespDiff;
imgStats.lateCycRespDiffxOriErr = oriGroups.lateCycRespDiffErr;
imgStats.lateCycRespTestxOri = oriGroups.lateCycTest;
imgStats.lateCycSITestxOri = oriGroups.lateCycSITestEaOri;
imgStats.lateCycSImatchedFRTestxOri = oriGroups.FRmatchedSITestEaOri;
% fprintf('Late Cyc Resp (Vis-Aud Trials) 0/45/90/135 pref = %s\n', ...
%     num2str(round(imgStats.lateCycRespDiffxOri,2,'significant')))
% fprintf('Late Cyc Resp Err (Vis-Aud Trials) 0/45/90/135 pref = %s\n', ...
%     num2str(round(imgStats.lateCycRespDiffxOriErr,2,'significant')))
fprintf('Late Cyc Resp AV Test x Ori, 0/45/90/135 pref, p=%s\n',num2str(round(...
    imgStats.lateCycRespTestxOri,2,'significant')))
fprintf('Late Cyc SI Test x Ori, 0/45/90/135 pref, a=%s, p=%s\n',...
    num2str(round(0.05/(nOri-1),2,'significant')),num2str(round(...
    imgStats.lateCycSITestxOri,2,'significant')))
fprintf('Late Cyc FR Matched SI Test x Ori, 0/45/90/135 pref, a=%s, p=%s\n',...
    num2str(round(0.05/(nOri-1),2,'significant')),num2str(round(...
    imgStats.lateCycSImatchedFRTestxOri,2,'significant')))
imgStats.lateCycSIxOriANOVA = oriGroups.lateCycSITest;
fprintf('SI x Ori 1-Way ANOVA, p=%s\n',num2str(round(...
    imgStats.lateCycSIxOriANOVA,2,'significant')))
imgStats.lateCycSIxOriMatchFRTest = oriGroups.matchSIxFRTest;
fprintf('SI x Ori Match FR, ANOVA p=%s\n',num2str(round(...
    imgStats.lateCycSIxOriMatchFRTest,2,'significant')))

%%
if strcmp(ds,'FSAV_attentionV1')
    minCellN_SIFRmatch = 36;
else
    minCellN_SIFRmatch = 4;
end
aurocGroups = struct;
aurocGroups.name = {'Dist.';'Tar';'n.d.'};
aurocGroups.cmp(1).name = 'first:target';
aurocGroups.cmp(2).name = 'late:target';
for icmp = 1:2
    aurocGroups.cmp(icmp).n = nan(1,3);
    aurocGroups.cmp(icmp).n4SI = nan(1,3);
    aurocGroups.cmp(icmp).targetTC = cell(1,3);
    aurocGroups.cmp(icmp).targetTCErr = cell(1,3); 
    aurocGroups.cmp(icmp).targetTC_distOnly = cell(1,3);
    aurocGroups.cmp(icmp).targetTCErr_distOnly = cell(1,3);    
    aurocGroups.cmp(icmp).firstTC = cell(1,3);
    aurocGroups.cmp(icmp).firstTCErr = cell(1,3);
    aurocGroups.cmp(icmp).lateCycSI = nan(1,3);
    aurocGroups.cmp(icmp).lateCycSIErr = nan(1,3);
    aurocGroups.cmp(icmp).lateCycAVauROC = nan(1,3);
    aurocGroups.cmp(icmp).lateCycAVauROCErr = nan(1,3);
%     aurocGroups.cmp(icmp).lateCycSI_temp = nan(1,2);
%     aurocGroups.cmp(icmp).lateCycSIErr_temp = nan(1,2);
    aurocGroups.cmp(icmp).lateCycSITestData = cell(1,3);
    aurocGroups.cmp(icmp).matchlateCycSITestData = cell(1,3);
    aurocGroups.cmp(icmp).matchLateCycSI = nan(1,3);
    aurocGroups.cmp(icmp).matchLateCycSIErr = nan(1,3);
    aurocGroups.cmp(icmp).matchLateCycResp = nan(1,3);
    aurocGroups.cmp(icmp).matchLateCycRespErr = nan(1,3);
end
for igrp = 1:3
    if igrp == 3
        ind_first = ~cellInfo.firstStimAuROCTest & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
        ind_late = ~cellInfo.lateStimAuROCTest & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    elseif igrp == 1
        ind_first = cellInfo.firstStimAuROCTest & cellInfo.firstStimAuROC < 0.5 & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
        ind_late = cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC < 0.5 & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    elseif igrp == 2
        ind_first = cellInfo.firstStimAuROCTest & cellInfo.firstStimAuROC > 0.5 & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
        ind_late = cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC > 0.5 & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    end
    aurocGroups.cmp(1).n(igrp) = sum(ind_first);
    aurocGroups.cmp(2).n(igrp) = sum(ind_late);
    aurocGroups.cmp(1).n4SI(igrp) = sum(ind_first &...
        (cellInfo.lateCycRespCells));
    aurocGroups.cmp(2).n4SI(igrp) = sum(ind_late &...
        (cellInfo.lateCycRespCells));
    aurocGroups.cmp(1).targetTC{igrp} = mean(targetAnalysis.tc{3,visualTrials}...
        (:,ind_first),2);
    aurocGroups.cmp(1).targetTCErr{igrp} = ste(targetAnalysis.tc{3,visualTrials}...
        (:,ind_first),2);
    aurocGroups.cmp(2).targetTC{igrp} = mean(targetAnalysis.tc{3,visualTrials}...
        (:,ind_late),2);
    aurocGroups.cmp(2).targetTCErr{igrp} = ste(targetAnalysis.tc{3,visualTrials}...
        (:,ind_late),2);
    aurocGroups.cmp(1).targetTC_distOnly{igrp} = mean(targetAnalysis.tc{3,visualTrials}...
        (:,ind_first & (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(1).targetTCErr_distOnly{igrp} = ste(targetAnalysis.tc{3,visualTrials}...
        (:,ind_first & (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).targetTC_distOnly{igrp} = mean(targetAnalysis.tc{3,visualTrials}...
        (:,ind_late & (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).targetTCErr_distOnly{igrp} = ste(targetAnalysis.tc{3,visualTrials}...
        (:,ind_late & (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(1).firstTC{igrp} = mean(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_first &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(1).firstTCErr{igrp} = ste(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_first &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).firstTC{igrp} = mean(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_late &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).firstTCErr{igrp} = ste(antiAnalysis.firstCycTC{visualTrials}...
        (:,ind_late &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(1).lateCycSI(igrp) = mean(antiAnalysis.lateCycSI(ind_first &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(1).lateCycSIErr(igrp) = ste(antiAnalysis.lateCycSI(ind_first &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).lateCycSI(igrp) = mean(antiAnalysis.lateCycSI(ind_late &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).lateCycSIErr(igrp) = ste(antiAnalysis.lateCycSI(ind_late &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).lateCycSITestData{igrp} = antiAnalysis.lateCycSI(ind_late &...
        (cellInfo.lateCycRespCells));
    aurocGroups.cmp(2).lateCycAVauROC(igrp) = mean(antiAnalysis.lateCycAVauROC(ind_late &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).lateCycAVauROCErr(igrp) = ste(antiAnalysis.lateCycAVauROC(ind_late &...
        (cellInfo.lateCycRespCells)),2);
    aurocGroups.cmp(2).lateCycAVauROCTestData{igrp} = antiAnalysis.lateCycAVauROC(ind_late &...
        (cellInfo.lateCycRespCells));
    
    lateCycResp_grp = lateCycRespAll(ind_late & cellInfo.lateCycRespCells);
    lateCycSI_grp = antiAnalysis.lateCycSI(ind_late & cellInfo.lateCycRespCells);
    [sortFR,sortFRInd] = sort(lateCycResp_grp,'descend');
    sortSI = lateCycSI_grp(sortFRInd);
    aurocGroups.cmp(2).matchLateCycSITestData{igrp} = sortSI(1:minCellN_SIFRmatch);
    aurocGroups.cmp(2).matchLateCycSI(igrp) = mean(sortSI(1:minCellN_SIFRmatch),2);
    aurocGroups.cmp(2).matchLateCycSIErr(igrp) = ste(sortSI(1:minCellN_SIFRmatch),2);
    aurocGroups.cmp(2).matchLateCycResp(igrp) = mean(sortFR(1:minCellN_SIFRmatch),2);
    aurocGroups.cmp(2).matchLateCycRespErr(igrp) = ste(sortFR(1:minCellN_SIFRmatch),2);
end
ind = (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
[~,sortInd] = sort(cellInfo.firstStimAuROC(ind));
tc = targetAnalysis.tc{3,visualTrials}(:,ind)';
aurocGroups.cmp(1).hm = flipud(tc(sortInd,:));
[~,sortInd] = sort(cellInfo.lateStimAuROC(ind));
tc = targetAnalysis.tc{3,visualTrials}(:,ind)';
aurocGroups.cmp(2).hm = flipud(tc(sortInd,:));

grpID = [];
for igrp = 1:3
    grpID = cat(2,grpID,ones(1,length(aurocGroups.cmp(2).lateCycSITestData{igrp})).*igrp);
end
[aurocGroups.cmp(2).lateCycSITest,~,stats] = anova1(cell2mat(...
    aurocGroups.cmp(2).lateCycSITestData),grpID,'off');
aurocGroups.cmp(2).lateCycSIPostHoc = multcompare(stats);
grpID = [];
for igrp = 1:3
    grpID = cat(2,grpID,ones(1,length(aurocGroups.cmp(2).matchLateCycSITestData{igrp})).*igrp);
end
[aurocGroups.cmp(2).matchLateCycSITest,~,stats] = anova1(cell2mat(...
    aurocGroups.cmp(2).matchLateCycSITestData),grpID,'off');
aurocGroups.cmp(2).lateCycSIMatchFRPostHoc = multcompare(stats);

setFigParams4Print('portrait')
figure
suptitle('Target Aligned Resp., All Task Resp. Cells')
colormap(brewermap([],'*RdBu'))
for icmp = 1:2
    subplot(2,2,icmp)
    imagesc(aurocGroups.cmp(1).hm)
    caxis(hmLim)
    colorbar
    figXAxis([],'Time (ms)',[],ttLabelFr_target,ttLabel_target)
    figYAxis([],'Cell # (auROC Sorted)',[])
    figAxForm
    title(aurocGroups.cmp(icmp).name)
%     subplot(3,2,icmp+2)
%     L = [];
%     for igrp = 1:3
%         y = aurocGroups.cmp(icmp).targetTC{igrp};
%         yerr = aurocGroups.cmp(icmp).targetTCErr{igrp};
%         hold on
%         h = shadedErrorBar_chooseColor(tt_targetTC,y,yerr,aurocColor{igrp});
%         L(igrp) = h.mainLine;
%     end
%     figXAxis([],'Time (ms)',[tt_targetTC(1) tt_targetTC(end)],...
%         ttLabel_target,ttLabel_target)
%     figYAxis([],'dF/F',targetTCLim)  
%     hline(0,'k:')
%     vline(respWinTT,'k--')
%     vline(preTargetStimLabel,'k:')
%     figAxForm
%     title(sprintf('%s, n=%s',aurocGroups.cmp(icmp).name,...
%         num2str(sum(aurocGroups.cmp(icmp).n))))
%     legend(L,strcat(aurocGroups.name,repmat(' n=',[3,1]),...
%         cellfun(@num2str,num2cell(aurocGroups.cmp(icmp).n)','unif',0)))

    subplot(2,2,icmp+2)
    for igrp = 1:3
        y = aurocGroups.cmp(icmp).firstTC{igrp}(tcStartFrame:end);
        yerr = aurocGroups.cmp(icmp).firstTCErr{igrp}(tcStartFrame:end);
        hold on
        h = shadedErrorBar_chooseColor(tt_cycTC,y,yerr,aurocColor{igrp});
    end
    figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
    figYAxis([],'dF/F',targetTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    figAxForm
    title(sprintf('%s, First Stim Resp.',aurocGroups.cmp(icmp).name))
end
print([fnout 'tuningAnalysis_aurocGroupsTC'],'-dpdf','-fillpage')  

setFigParams4Print('landscape')
figure
suptitle('Distractor Resp. Cells')
for icmp = 1:2
    subplot(3,2,icmp)
    L = [];
    for igrp = 1:3
        y = aurocGroups.cmp(icmp).targetTC_distOnly{igrp};
        yerr = aurocGroups.cmp(icmp).targetTCErr_distOnly{igrp};
        hold on
        h = shadedErrorBar_chooseColor(tt_targetTC,y,yerr,aurocColor{igrp});
        L(igrp) = h.mainLine;
    end
    figXAxis([],'Time (ms)',[-700 350],...
        ttLabel_target,ttLabel_target)
    figYAxis([],'dF/F',targetTCLim)  
    hline(0,'k:')
    vline(respWinTT,'k--')
    vline(preTargetStimLabel,'k:')
    figAxForm
    title(sprintf('%s, n=%s',aurocGroups.cmp(icmp).name,...
        num2str(sum(aurocGroups.cmp(icmp).n4SI))))
    legend(L,strcat(aurocGroups.name,repmat(' n=',[3,1]),...
        cellfun(@num2str,num2cell(aurocGroups.cmp(icmp).n4SI)','unif',0)))
    
    subplot(3,2,icmp+2)
%     bar(aurocGroups.cmp(icmp).lateCycSI,'EdgeColor','none','BarWidth',0.5);
    hold on
    errorbar(1:3,aurocGroups.cmp(icmp).lateCycSI,aurocGroups.cmp(icmp).lateCycSIErr,'.')
    figXAxis([],'Pref. Task Stim. (by auROC)',[0 4],1:3,aurocGroups.name)
    figYAxis([],'Selectivity',siOriLim)
    figAxForm
    hline(0,'k:')
    for igrp = 1:3
        text(igrp,aurocGroups.cmp(icmp).lateCycSI(igrp)+.2,num2str(aurocGroups.cmp(icmp).n4SI(igrp)));
    end
    title(aurocGroups.cmp(icmp).name)
    
    if icmp == 2
        subplot(3,2,icmp+4)
        hold on
        errorbar(1:3,aurocGroups.cmp(icmp).lateCycAVauROC,...
            aurocGroups.cmp(icmp).lateCycAVauROCErr,'.')
        figXAxis([],'Pref. Task Stim. (by auROC)',[0 4],1:3,aurocGroups.name)
        figYAxis([],'auROC_V_A',[0.45 0.55])
        figAxForm
        hline(0,'k:')
%         for igrp = 1:3
%             text(igrp,aurocGroups.cmp(icmp).lateCycSI(igrp)+.2,num2str(aurocGroups.cmp(icmp).n4SI(igrp)));
%         end
        title(aurocGroups.cmp(icmp).name)
    end
end
print([fnout 'tuningAnalysis_aurocGroupsSI'],'-dpdf','-fillpage')

%% task tuning binned by SI (Figure 3)
ind = cellInfo.lateCycRespCells & antiAnalysis.lateCycAV95CITest';
frTTInBin = nan(nSIBins,3);
frSIInTTBin = nan(3,nSIBins);
for ibin = 1:nSIBins
    binInd = siBinID' == ibin & ind;
    frTTInBin(ibin,:) = histcounts(aurocBinID(binInd),3,...
        'Normalization','probability');
%     frTTInBin(ibin,:) = histcounts(aurocBinID(binInd),3);
end
for itt = 1:3
    binInd = ind & aurocBinID' == itt;
    frSIInTTBin(itt,:) = histcounts(siBinID(binInd),nSIBins,...
        'Normalization','probability');
end
binnedHT = cell(3,nSIBins);
binnedET = cell(3,nSIBins);
binnedT = cell(3,nSIBins);
binnedD = cell(3,nSIBins);
binnedSI = cell(3,nSIBins);
for ibin = 1:nSIBins
    for itt = 1:3
        binInd = ind & siBinID' == ibin & aurocBinID' == itt;
        binnedHT{itt,ibin} = Tr_vis{1}(binInd);
        binnedET{itt,ibin} = Tr_vis{2}(binInd);
        binnedT{itt,ibin} = Tr_vis{3}(binInd);
        binnedD{itt,ibin} = Dr{3}(binInd);
        binnedSI{itt,ibin} = SI(binInd);
    end
end

matName = {'Stim';'SI'};
ttVisSIBinTest = cell(1,3);
for itt = 1:3
    d = cat(1,binnedD(itt,:),binnedT(itt,:));
%     d = cat(1,binnedD(itt,:),binnedHT(itt,:));
    ttVisSIBinTest{itt} = twoWayAnovaUnmatched(d,matName);
end
n = cellfun(@length,binnedD);


for itt = 1:3
    fprintf('%s pref - SI x D/T Two-Way ANOVA, p: %s=%s, %s=%s int=%s\n',...
        aurocGroups.name{itt},...
        ttVisSIBinTest{itt}.varTestName{1},num2str(round(ttVisSIBinTest{itt}.p(1),2,'significant')),...
        ttVisSIBinTest{itt}.varTestName{2},num2str(round(ttVisSIBinTest{itt}.p(2),2,'significant')),...
        num2str(round(ttVisSIBinTest{itt}.p(3),2,'significant')))
    if any(ttVisSIBinTest{itt}.p < 0.05)
        for ip = 1:2
            disp(ttVisSIBinTest{itt}.posthoc{ip})
        end
    end
end

% Supplemental Fig 6: distractor response distribution
distRespFig = figure;
subplot 121
L = [];
cAll = get(gca,'colororder');
dFFDist = cell(1,3);
for itt = 1:3
    d = cell2mat(binnedD(itt,:));
    dFFDist{itt} = d;
    hold on
    h = cdfplot(d);
    h.Color = cAll(itt,:);  
    h2 = vline(mean(d),'-');
    h2.Color = h.Color;
    L(itt) = h;
end
legend(L,{'D','T','NP'},...
    'location','northeastoutside')
figXAxis([],'Dist. Resp (dF/F)',[0 0.04])
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
[p,~,stats] = anova1(Dr{3}(ind),aurocBinID(ind),'off');
ph = multcompare(stats,'display','off');
if p < 0.05
    fprintf('Vis Task Tuning Dist Resp Distributions:\n')
    disp(ph(:,[1:2,6]))
end
title(sprintf('Vis task tuning,p=%s,',num2str(round(p,2,'significant'))))

figure
suptitle('Visual Trial Task Tuning, Tuned Late Resp Cells')
subplot 321
bar(frTTInBin,'stacked')
figXAxis([],'',[0 4],1:3,{'-SI','0','+SI'})
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
legend({'D','T','NP'},...
    'location','northeastoutside')
subplot 322
bar(frSIInTTBin,'stacked')
figXAxis([],'',[0 4],1:3,{'D','T','NP'})
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
legend({'-SI','0','+SI'},...
    'location','northeastoutside')
for itt = 1:3
    subplot(3,2,itt+2)
    x = cellfun(@mean,binnedSI(itt,:));
    xerr = cellfun(@(x) ste(x,2),binnedSI(itt,:));
    y = cellfun(@mean,binnedD(itt,:));
    yerr = cellfun(@(x) ste(x,2),binnedD(itt,:));    
    hold on
    errorbar(x,y,yerr,yerr,xerr,xerr,'.-')      
%     y = cellfun(@mean,binnedHT(itt,:));
%     yerr = cellfun(@(x) ste(x,2),binnedHT(itt,:));    
%     errorbar(x,y,yerr,yerr,xerr,xerr,'.-')      
%     y = cellfun(@mean,binnedET(itt,:));
%     yerr = cellfun(@(x) ste(x,2),binnedET(itt,:));    
    y = cellfun(@mean,binnedT(itt,:));
    yerr = cellfun(@(x) ste(x,2),binnedT(itt,:));        hold on
    errorbar(x,y,yerr,yerr,xerr,xerr,'.-')   
    figXAxis([],'SI (mean in bin)',siBinLim)
    figYAxis([],'Stim Resp (dF/F)',stimRespLim)
    figAxForm   
%     L=legend({'D_a_l_l','HT','ET'},'location','northeastoutside');
    L=legend({'D_a_l_l','T_a_l_l'},'location','northeastoutside');
end
subplot 323
title(sprintf('Pref. Dist., n=%s',num2str(sum(n(1,:)))))
subplot 324
title(sprintf('Pref. Target., n=%s',num2str(sum(n(2,:)))))
subplot 325
title(sprintf('No Pref., n=%s',num2str(sum(n(3,:)))))
print([fnout 'tuningAnalysis_aurocGroups_SIBinned'],'-dpdf','-fillpage')

frTTInBin = nan(nSIBins,3);
frSIInTTBin = nan(3,nSIBins);
for ibin = 1:nSIBins
    binInd = siBinID' == ibin & ind;
    frTTInBin(ibin,:) = histcounts(aurocBinID_aud(binInd),3,...
        'Normalization','probability');
%     frTTInBin(ibin,:) = histcounts(aurocBinID(binInd),3);
end
for itt = 1:3
    binInd = ind & aurocBinID_aud' == itt;
    frSIInTTBin(itt,:) = histcounts(siBinID(binInd),nSIBins,...
        'Normalization','probability');
end
binnedT = cell(3,nSIBins);
binnedD = cell(3,nSIBins);
binnedSI = cell(3,nSIBins);
for ibin = 1:nSIBins
    for itt = 1:3
        binInd = ind & siBinID' == ibin & aurocBinID_aud' == itt;
        binnedT{itt,ibin} = Tr_aud{3}(binInd);
        binnedD{itt,ibin} = Dr{3}(binInd);
        binnedSI{itt,ibin} = SI(binInd);
    end
end

matName = {'Stim';'SI'};
ttAudSIBinTest = cell(1,3);
for itt = 1:3
    d = cat(1,binnedD(itt,:),binnedT(itt,:));
%     d = cat(1,binnedD(itt,:),binnedHT(itt,:));
    ttAudSIBinTest{itt} = twoWayAnovaUnmatched(d,matName);
end
n = cellfun(@length,binnedD);

figure
suptitle('Auditory Trial Task Tuning')
subplot 321
bar(frTTInBin,'stacked')
figXAxis([],'',[0 4],1:3,{'-SI','0','+SI'})
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
legend({'D','T','NP'},...
    'location','northeastoutside')
subplot 322
bar(frSIInTTBin,'stacked')
figXAxis([],'',[0 4],1:3,{'D','T','NP'})
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
legend({'-SI','0','+SI'},...
    'location','northeastoutside')
for itt = 1:3
    subplot(3,2,itt+2)
    x = cellfun(@mean,binnedSI(itt,:));
    xerr = cellfun(@(x) ste(x,2),binnedSI(itt,:));
    y = cellfun(@mean,binnedD(itt,:));
    yerr = cellfun(@(x) ste(x,2),binnedD(itt,:));    
    hold on
    errorbar(x,y,yerr,yerr,xerr,xerr,'.-')      
    y = cellfun(@mean,binnedT(itt,:));
    yerr = cellfun(@(x) ste(x,2),binnedT(itt,:));    
    hold on
    errorbar(x,y,yerr,yerr,xerr,xerr,'.-')   
    figXAxis([],'SI (mean in bin)',siBinLim)
    figYAxis([],'Stim Resp (dF/F)',stimRespLim)
    figAxForm   
    L=legend({'D_a_l_l','T_a_l_l'},'location','northeastoutside');
end
subplot 323
title(sprintf('Pref. Dist., n=%s',num2str(sum(n(1,:)))))
subplot 324
title(sprintf('Pref. Target, n=%s',num2str(sum(n(2,:)))))
subplot 325
title(sprintf('No Pref., n=%s',num2str(sum(n(3,:)))))
print([fnout 'tuningAnalysis_audAurocGroups_SIBinned'],'-dpdf','-fillpage')

figure(distRespFig);
subplot 122
L = [];
cAll = get(gca,'colororder');
dFFDist = cell(1,3);
for itt = 1:3
    d = cell2mat(binnedD(itt,:));
    dFFDist{itt} = d;
    hold on
    h = cdfplot(d);
    h.Color = cAll(itt,:);  
    h2 = vline(mean(d),'-');
    h2.Color = h.Color;
    L(itt) = h;
end
legend(L,{'D','T','NP'},...
    'location','northeastoutside')
figXAxis([],'Dist. Resp (dF/F)',[0 0.04])
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
[p,~,stats] = anova1(Dr{3}(ind),aurocBinID_aud(ind),'off');
ph = multcompare(stats,'display','off');
if p < 0.05
    fprintf('Aud Task Tuning Dist Resp Distributions:\n')
    disp(ph(:,[1:2,6]))
end
title(sprintf('Aud task tuning,p=%s,',num2str(round(p,2,'significant'))))
print([fnout 'tuningAnalysis_distRespDistrib_aurocGroups'],'-dpdf','-fillpage')
%% auditory auroc groups comparison

if strcmp(ds,'FSAV_attentionV1')%%
minCellN_SIFRmatch = 61;
maxFR_lowgrp = 0.0264;
audAurocGroups = struct;
for igrp = 1:3
    if igrp == 3
        ind_late = ~cellInfo.audLateStimAuROCTest & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    elseif igrp == 1
        ind_late = cellInfo.audLateStimAuROCTest & cellInfo.audLateStimAuROC < 0.5 & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    elseif igrp == 2
        ind_late = cellInfo.audLateStimAuROCTest & cellInfo.audLateStimAuROC > 0.5 & ...
            (cellInfo.lateCycRespCells | cellInfo.targetRespCells);
    end
    audAurocGroups.n(igrp) = sum(ind_late);
    audAurocGroups.n4SI(igrp) = sum(ind_late & cellInfo.lateCycRespCells);
%     audAurocGroups.targetTC_distOnly{igrp} = mean(targetAnalysis.tc{3,auditoryTrials}...
%         (:,ind_late & (cellInfo.lateCycRespCells)),2);
%     audAurocGroups.targetTCErr_distOnly{igrp} = ste(targetAnalysis.tc{3,auditoryTrials}...
%         (:,ind_late & (cellInfo.lateCycRespCells)),2);
    audAurocGroups.lateCycSI(igrp) = mean(antiAnalysis.lateCycSI(ind_late &...
        (cellInfo.lateCycRespCells)),2);
    audAurocGroups.lateCycSIErr(igrp) = ste(antiAnalysis.lateCycSI(ind_late &...
        (cellInfo.lateCycRespCells)),2);
    audAurocGroups.lateCycSITestData{igrp} = antiAnalysis.lateCycSI(ind_late &...
        (cellInfo.lateCycRespCells));
    audAurocGroups.targetTC_distOnly{igrp} = mean(targetAnalysis.tc{3,auditoryTrials}...
        (:,ind_late & (cellInfo.lateCycRespCells)),2);
    audAurocGroups.targetTCErr_distOnly{igrp} = ste(targetAnalysis.tc{3,auditoryTrials}...
        (:,ind_late & (cellInfo.lateCycRespCells)),2);
    lateCycResp_grp = lateCycRespAll(ind_late & cellInfo.lateCycRespCells);
    lateCycSI_grp = antiAnalysis.lateCycSI(ind_late & cellInfo.lateCycRespCells);
    [sortFR,sortFRInd] = sort(lateCycResp_grp,'descend');
    sortSI = lateCycSI_grp(sortFRInd); 
    if igrp == 1
%         ind1 = find(sortFR <= maxFR_lowgrp,1);
        ind1 = 10;
        matchInd = ind1:(minCellN_SIFRmatch+ind1-1);
    elseif igrp == 2
        matchInd = 1:minCellN_SIFRmatch;
    elseif igrp == 3
%         ind1 = find(sortFR <= maxFR_lowgrp,1);
        ind1 = 10;
        matchInd = ind1:2:2*(minCellN_SIFRmatch+ind1-1);
    end
    audAurocGroups.matchLateCycSITestData{igrp} = sortSI(matchInd);
    audAurocGroups.matchLateCycSI(igrp) = mean(sortSI(matchInd),2);
    audAurocGroups.matchLateCycSIErr(igrp) = ste(sortSI(matchInd),2);
    audAurocGroups.matchLateCycResp(igrp) = mean(sortFR(matchInd),2);
    audAurocGroups.matchLateCycRespErr(igrp) = ste(sortFR(matchInd),2);
end

grpID = [];
for igrp = 1:3
    grpID = cat(2,grpID,ones(1,length(audAurocGroups.lateCycSITestData{igrp})).*igrp);
end
[audAurocGroups.lateCycSITest,~,stats] = anova1(cell2mat(...
    audAurocGroups.lateCycSITestData),grpID,'off');
audAurocGroups.lateCycSIPostHoc = multcompare(stats,[],'off');
grpID = [];
for igrp = 1:3
    grpID = cat(2,grpID,ones(1,length(audAurocGroups.matchLateCycSITestData{igrp})).*igrp);
end
[audAurocGroups.matchLateCycSITest,~,stats] = anova1(cell2mat(...
    audAurocGroups.matchLateCycSITestData),grpID,'off');
audAurocGroups.matchLateCycSIPostHoc = multcompare(stats);

figure
subplot 121
L = [];
for igrp = 1:3
    y = audAurocGroups.targetTC_distOnly{igrp};
    yerr = audAurocGroups.targetTCErr_distOnly{igrp};
    hold on
    h = shadedErrorBar_chooseColor(tt_targetTC,y,yerr,aurocColor{igrp});
    L(igrp) = h.mainLine;
end
figXAxis([],'Time (ms)',[-700 350],...
        ttLabel_target,ttLabel_target)
figYAxis([],'dF/F',targetTCLim)  
hline(0,'k:')
vline(respWinTT,'k--')
vline(preTargetStimLabel,'k:')
figAxForm
title(sprintf('%s, n=%s',aurocGroups.cmp(icmp).name,...
    num2str(sum(audAurocGroups.n4SI))))
legend(L,strcat(aurocGroups.name,repmat(' n=',[3,1]),...
    cellfun(@num2str,num2cell(audAurocGroups.n4SI)','unif',0)))
subplot 122
% bar(audAurocGroups.lateCycSI,'EdgeColor','none','BarWidth',0.5);hold on
errorbar(1:3,audAurocGroups.lateCycSI,audAurocGroups.lateCycSIErr,'.')
figXAxis([],'Pref. Aud. Task Stim. (by auROC)',[0 4],1:3,aurocGroups.name)
figYAxis([],'Selectivity',siOriLim)
figAxForm
hline(0,'k:')
for igrp = 1:3
    text(igrp,audAurocGroups.lateCycSI(igrp)+.2,num2str(audAurocGroups.n4SI(igrp)));
end
print([fnout 'tuningAnalysis_audAurocGroupsSI'],'-dpdf')
%% Vis and Aud Selectivity matched by FR
figure
suptitle('Late Stim auROC groups, SI matched by late stim FR on vis trials')
subplot 221
% bar(aurocGroups.cmp(2).matchLateCycSI,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:3,aurocGroups.cmp(2).matchLateCycSI,aurocGroups.cmp(2).matchLateCycSIErr,'.')
figXAxis([],'Pref. Task Stim. (by auROC)',[0 4],1:3,aurocGroups.name)
figYAxis([],'Selectivity',siOriLim)
figAxForm
hline(0,'k:')
title(sprintf('Visual auROC, anova=%s',num2str(aurocGroups.cmp(2).matchLateCycSITest)))
subplot 222
% bar(audAurocGroups.matchLateCycSI,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:3,audAurocGroups.matchLateCycSI,audAurocGroups.matchLateCycSIErr,'.')
figXAxis([],'Pref. Task Stim. (by auROC)',[0 4],1:3,aurocGroups.name)
figYAxis([],'Selectivity',siOriLim)
figAxForm
hline(0,'k:')
title(sprintf('Auditory auROC, anova=%s',num2str(audAurocGroups.matchLateCycSITest)))
subplot 223
% bar(aurocGroups.cmp(2).matchLateCycResp,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:3,aurocGroups.cmp(2).matchLateCycResp,...
    aurocGroups.cmp(2).matchLateCycRespErr,'.')
figXAxis([],'Pref. Task Stim. (by auROC)',[0 4],1:3,aurocGroups.name)
figYAxis([],'dF/F',[0 0.02])
figAxForm
subplot 224
% bar(audAurocGroups.matchLateCycResp,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:3,audAurocGroups.matchLateCycResp,...
    audAurocGroups.matchLateCycRespErr,'.')
figXAxis([],'Pref. Task Stim. (by auROC)',[0 4],1:3,aurocGroups.name)
figYAxis([],'dF/F',[0 0.02])
figAxForm

print([fnout 'tuningAnalysis_aurocGroupsSImatchFR'],'-dpdf','-fillpage')

end
%% stats for task tuning groups
imgStats.nCells.aurocGroups = aurocGroups.cmp(2).n4SI;
fprintf('n D/T/NP pref = %s\n', num2str(imgStats.nCells.aurocGroups))

imgStats.lateCycSIxAGTest = aurocGroups.cmp(2).lateCycSITest;
imgStats.lateCycSIxAGPostHocTest = aurocGroups.cmp(2).lateCycSIPostHoc(:,[1,2,end]);
fprintf('Late Cyc SI x Task Tuning, ANOVA p=%s\n',...
    num2str(imgStats.lateCycSIxAGTest))
disp(imgStats.lateCycSIxAGPostHocTest)

imgStats.lateCycFRMatchSIxAGTest = aurocGroups.cmp(2).matchLateCycSITest;
imgStats.lateCycFRMatchSIxAGPostHocTest = aurocGroups.cmp(2).lateCycSIMatchFRPostHoc(:,[1,2,end]);
fprintf('Late Cyc Matched FR SI x Task Tuning, ANOVA p=%s\n',...
    num2str(imgStats.lateCycFRMatchSIxAGTest))
disp(imgStats.lateCycFRMatchSIxAGPostHocTest)

%% SI logarithmically fit to FR (Supplemental Fig. 7)

rng(0)
nBoot = 1000;
aurocGroups2 = struct;
aurocGroups2.name = {'T or D';'NP'};
corrBoot = nan(nBoot,2);
rsqBoot = nan(nBoot,2);
for iboot = 1:nBoot    
    ind = find(cellInfo.lateCycRespCells & cellInfo.lateStimAuROCTest);
    sampInd = randsample(ind,200,1);
    x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,sampInd),1);
    x = log(x);
    y = antiAnalysis.lateCycSI(sampInd);
    mdl = fitlm(x,y);
    rsqBoot(iboot,1) = mdl.Rsquared.Ordinary;
    corrmat = corrcoef(x,y);
    corrBoot(iboot,1) = corrmat(1,2);
    
    ind = find(cellInfo.lateCycRespCells & ~cellInfo.lateStimAuROCTest);
    sampInd = randsample(ind,200,1);
    x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,sampInd),1);
    x = log(x);
    y = antiAnalysis.lateCycSI(sampInd);
    mdl = fitlm(x,y);
    rsqBoot(iboot,2) = mdl.Rsquared.Ordinary;
    corrmat = corrcoef(x,y);
    corrBoot(iboot,2) = corrmat(1,2);
end
ciInd = [round(nBoot*0.05), round(nBoot*0.95)];
aurocGroups2.siFRCorr = mean(corrBoot,1);
corrBootSort = [sort(corrBoot(:,1)), sort(corrBoot(:,2))];
aurocGroups2.siFRCorr95ci = nan(2,2);
aurocGroups2.siFRCorr95ci(:,1) = corrBootSort(ciInd,1);
aurocGroups2.siFRCorr95ci(:,2) = corrBootSort(ciInd,2);
ind = cellInfo.lateCycRespCells & cellInfo.lateStimAuROCTest;
% x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1);
% y = antiAnalysis.lateCycSI(ind);
% corrmat = corrcoef(x,y);
% c_TD = corrmat(1,2);
n_TD = sum(ind);
ind = cellInfo.lateCycRespCells & ~cellInfo.lateStimAuROCTest;
% x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1);
% y = antiAnalysis.lateCycSI(ind);
% corrmat = corrcoef(x,y);
% c_NP = corrmat(1,2);
n_NP = sum(ind);
% aurocGroups2.siFRCorr = [c_TD,c_NP];
aurocGroups2.siFRCorrTest = compare_correlation_coefficients(...
    aurocGroups2.siFRCorr(1),aurocGroups2.siFRCorr(2),n_TD,n_NP);
aurocGroups2.siFRRsq = mean(rsqBoot,1);
rsqBootSort = [sort(rsqBoot(:,1)), sort(rsqBoot(:,2))];
aurocGroups2.siFRRsq95ci = nan(2,2);
aurocGroups2.siFRRsq95ci(:,1) = rsqBootSort(ciInd,1);
aurocGroups2.siFRRsq95ci(:,2) = rsqBootSort(ciInd,2);

respBinEdges = [-1, -0.002:0.002:0.006, 1];
figure
subplot 231
% ind = (cellInfo.lateCycRespCells) & ...
%     mean(antiAnalysis.lateCycTC{visualTrials}(respwin,:),1)' < 0.032;
ind = (cellInfo.lateCycRespCells);
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1);
x = log(x);
y = antiAnalysis.lateCycSI(ind);
plot(x,y,'.')
mdl = fitlm(x,y);
yfit = predict(mdl,sort(x)');
hold on
plot(sort(x),yfit,'-')
[corrmat,pmat] = corrcoef(x,y);
figXAxis([],'Late Vis Stim. Resp.(log(dF/F))',[-8, -2])
ax = gca;
% ax.XScale = 'log';
figYAxis([],'Selectivity',siLim)
figAxForm
title({'All Late Stim Resp Neurons';sprintf('Linear Fit: R=%s,Rsq=%s,p=%s',...
    num2str(round(corrmat(1,2),2,'significant')),...
    num2str(round(mdl.Rsquared.Ordinary,2,'significant')),...
    num2str(round(pmat(1,2),2,'significant')))})

subplot 233
ind = cellInfo.lateCycRespCells & ~cellInfo.lateStimAuROCTest;
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1);
x = log(x);
y = antiAnalysis.lateCycSI(ind);
plot(x,y,'.')
mdl = fitlm(x,y);
yfit = predict(mdl,sort(x)');
hold on
plot(sort(x),yfit,'-')
[corrmat,pmat] = corrcoef(x,y);
figXAxis([],'Late Vis Stim. Resp.(log(dF/F))',[-8, -2])
ax = gca;
% ax.XScale = 'log';
figYAxis([],'Selectivity',siLim)
figAxForm
title({'No Pref. Neurons'; sprintf('Linear Fit: R=%s,Rsq=%s,p=%s',...
    num2str(round(corrmat(1,2),2,'significant')),...
    num2str(round(mdl.Rsquared.Ordinary,2,'significant')),...
    num2str(round(pmat(1,2),2,'significant')))})

subplot 232
ind = cellInfo.lateCycRespCells & cellInfo.lateStimAuROCTest;
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1);
x = log(x);
y = antiAnalysis.lateCycSI(ind);
mdl = fitlm(x,y);
yfit = predict(mdl,sort(x)');
hold on
plot(sort(x),yfit,'-')
[corrmat,pmat] = corrcoef(x,y);
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
    ind & cellInfo.lateStimAuROC < 0.5),1);
x = log(x);
y = antiAnalysis.lateCycSI(ind & cellInfo.lateStimAuROC < 0.5);
plot(x,y,'.')
x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
    ind & cellInfo.lateStimAuROC > 0.5),1);
x = log(x);
y = antiAnalysis.lateCycSI(ind & cellInfo.lateStimAuROC > 0.5);
plot(x,y,'.')
figXAxis([],'Late Vis Stim. Resp.(log(dF/F))',[-8, -2])
ax = gca;
% ax.XScale = 'log';
figYAxis([],'Selectivity',siLim)
figAxForm
title({'Dist. or Tar. Pref. Neurons'; sprintf('Linear Fit: R=%s,Rsq=%s,p=%s',...
    num2str(round(corrmat(1,2),2,'significant')),...
    num2str(round(mdl.Rsquared.Ordinary,2,'significant')),...
    num2str(round(pmat(1,2),2,'significant')))})
legend({'fit','Pref. Dist.','Pref. Tar.'},'location','southeast')

subplot 234
bar(aurocGroups2.siFRCorr,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:2,aurocGroups2.siFRCorr,aurocGroups2.siFRCorr - aurocGroups2.siFRCorr95ci(1,:),...
    aurocGroups2.siFRCorr95ci(2,:) - aurocGroups2.siFRCorr,'.');
figXAxis([],'',[0 3],1:2,aurocGroups2.name)
figYAxis([],'Bootstrapped Correlation (R)',[])
figAxForm
title(sprintf('p=%s, Fishers r to z trnsfm', ...
    num2str(round(aurocGroups2.siFRCorrTest,2,'significant'))))

subplot 235
bar(aurocGroups2.siFRRsq,'EdgeColor','none','BarWidth',0.5);
hold on
errorbar(1:2,aurocGroups2.siFRRsq,aurocGroups2.siFRRsq - aurocGroups2.siFRRsq95ci(1,:),...
    aurocGroups2.siFRRsq95ci(2,:) - aurocGroups2.siFRRsq,'.');
figXAxis([],'',[0 3],1:2,aurocGroups2.name)
figYAxis([],'Bootstrapped Rsq',[])
figAxForm

print([fnout 'siXVisTrialsdFFCorr_lateRespCells'],'-dpdf','-fillpage')
%% suppressed cells analysis (Supplemental Figure 5)
setFigParams4Print('portrait')
figure
subplot 311
ind = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells);
for iav = 1:2
    y = mean(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(ind) & ...
        cellInfo.isShortCycExpt),2);
    yerr = ste(antiAnalysis.longTC{iav}...
        ((tcStartFrame:end),(ind) & cellInfo.isShortCycExpt),2);
    hold on
    shadedErrorBar_chooseColor(tt_longTC,y,yerr,cueColor{iav});
end
figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
figYAxis([],'dF/F',suppTCLim)  
vline(lateWinTT,'k--')
hline(0,'k:')
figAxForm([],0)
title(sprintf('Suppressed Only Cells (%s/%s)',...
    num2str(sum((...
    ind)...
    & cellInfo.isShortCycExpt)),...
    num2str(sum(cellInfo.isShortCycExpt))))

subplot 323
x = mean(antiAnalysis.longTC{visualTrials}(lateWinFr,...
    (ind)),1);
xerr = ste(x,2);
y = mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,...
    (ind)),1);
yerr = ste(y,2);
plot(x,y,'.')
hold on
errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.')
plot(suppScatLim_win,suppScatLim_win,'k--')
plot(suppScatLim_win,[0 0],'k:')
plot([0 0],suppScatLim_win,'k:')
[~,p] = ttest(x,y);
figXAxis([],'Visual (dF/F)',suppScatLim_win)
figYAxis([],'Auditory (dF/F)',suppScatLim_win)
figAxForm
title(sprintf('Late Window, All Resp. Cells (%s/%s), p = %s',...
    num2str(sum(ind)),...
    num2str(length(cellInfo.firstRespCells)),...
    num2str(round(p,2,'significant'))))

subplot 324
y = antiAnalysis.lateWinSI(cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells);
[~,p]=ttest(y);
h = cdfplot(y);
hold on;
vline(mean(y),'k-')
figXAxis([],'Selectivity Index',siLim)
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
title(sprintf('Late Cyc. Resp. Cells, mean = %s+/-%s, p=%s',...
    num2str(round(mean(y),2,'significant')),...
    num2str(round(ste(y,2),2,'significant')),...
    num2str(round(p,2,'significant'))))

subplot 325
ind1 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC > 0.5;
ind2 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    cellInfo.lateStimAuROCTest & cellInfo.lateStimAuROC < 0.5;
ind3 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    ~cellInfo.lateStimAuROCTest;

y = antiAnalysis.lateWinSI(ind1);
yerr = ste(y,2);
errorbar(1,mean(y),yerr,'.')
[~,p]=ttest(y);
text(1,0.5,num2str(round(p,2,'significant')))
hold on
y = antiAnalysis.lateWinSI(ind2);
yerr = ste(y,2);
errorbar(2,mean(y),yerr,'.')
[~,p]=ttest(y);
text(2,0.5,num2str(round(p,2,'significant')))
y = antiAnalysis.lateWinSI(ind3);
yerr = ste(y,2);
errorbar(3,mean(y),yerr,'.')
[~,p]=ttest(y);
text(3,0.5,num2str(round(p,2,'significant')))
figXAxis([],'Vis. Task Stim. Pref.',[0 4],1:3,{'T','D','NP'})
figYAxis([],'Selectivity (Late Win)',[-1 0.5])
figAxForm
hline(0,'k:')
legend({num2str(sum(ind1)),num2str(sum(ind2)),num2str(sum(ind3))})

subplot 326
ind1 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    cellInfo.audLateStimAuROCTest & cellInfo.audLateStimAuROC > 0.5;
ind2 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    cellInfo.audLateStimAuROCTest & cellInfo.audLateStimAuROC < 0.5;
ind3 = cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells) & ...
    ~cellInfo.audLateStimAuROCTest;

y = antiAnalysis.lateWinSI(ind1);
yerr = ste(y,2);
errorbar(1,mean(y),yerr,'.')
[~,p]=ttest(y);
text(1,0.5,num2str(round(p,2,'significant')))
hold on
y = antiAnalysis.lateWinSI(ind2);
yerr = ste(y,2);
errorbar(2,mean(y),yerr,'.')
[~,p]=ttest(y);
text(2,0.5,num2str(round(p,2,'significant')))
y = antiAnalysis.lateWinSI(ind3);
yerr = ste(y,2);
errorbar(3,mean(y),yerr,'.')
[~,p]=ttest(y);
text(3,0.5,num2str(round(p,2,'significant')))
figXAxis([],'Aud. Task Stim. Pref.',[0 4],1:3,{'T','D','NP'})
figYAxis([],'Selectivity (Late Win)',[-1 0.5])
figAxForm
hline(0,'k:')
legend({num2str(sum(ind1)),num2str(sum(ind2)),num2str(sum(ind3))})

print([fnout 'tcLongTrialsAndLateCycWithQuant_SuppCells'],'-dpdf','-fillpage')
%% Supplemental Figure 6.
ind = antiAnalysis.lateCycAV95CITest == 1;
pctLateAVMod_UD = ([sum(ind' & antiAnalysis.lateCycSI' <= 0 & cellInfo.lateCycRespCells),...
    sum(ind' & antiAnalysis.lateCycSI' > 0 & cellInfo.lateCycRespCells)])...
    ./sum(cellInfo.lateCycRespCells);
pctSuppAVMod_UD = ([sum(ind' & antiAnalysis.lateWinSI' <= 0 &  cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells)),...
    sum(ind' & antiAnalysis.lateWinSI' > 0 & cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells)])...
    ./sum(cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells);

figure
subplot 211
bar(cat(1,pctLateAVMod_UD,pctSuppAVMod_UD),'stacked');
figXAxis([],'Cell Type',[0 3],1:2,{'Late Resp.';'Late Supp.'})
figYAxis([],'Fraction V-A Modulated Cells',[0 1])
figAxForm
legend({'V<A';'V>A'},'location','northeastoutside')
title('Significance tested by 95% CI of bootstrap V-A S.I.')

ind = antiAnalysis.lateCycAVShuffTest == 1;
pctLateAVMod_UD = ([sum(ind' & antiAnalysis.lateCycSI' <= 0 & cellInfo.lateCycRespCells),...
    sum(ind' & antiAnalysis.lateCycSI' > 0 & cellInfo.lateCycRespCells)])...
    ./sum(cellInfo.lateCycRespCells);
pctSuppAVMod_UD = ([sum(ind' & antiAnalysis.lateWinSI' <= 0 &  cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells)),...
    sum(ind' & antiAnalysis.lateWinSI' > 0 & cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells)])...
    ./sum(cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells);

subplot 212
bar(cat(1,pctLateAVMod_UD,pctSuppAVMod_UD),'stacked');
figXAxis([],'Cell Type',[0 3],1:2,{'Late Resp.';'Late Supp.'})
figYAxis([],'Fraction V-A Modulated Cells',[0 1])
figAxForm
legend({'V<A';'V>A'},'location','northeastoutside')
title('Significant SI - Shuffled Trial ID')

print([fnout 'fractionVAModCells'],'-dpdf')

%%
avName = {'Vis.','Aud'};
ind = cellInfo.lateSuppCells & ~cellInfo.lateCycRespCells;
for iav = 1:2
    imgStats.av(iav).lateWinSupp = mean(mean(antiAnalysis.longTC{iav}(lateWinFr,...
        (ind)),1));
    imgStats.av(iav).lateWinSuppErr = ste(mean(antiAnalysis.longTC{iav}(lateWinFr,...
        (ind)),1),2);
    fprintf('%s Suppr. Late Win mean/err: %s/%s\n',avName{iav},...
        num2str(round(imgStats.av(iav).lateWinSupp,2,'significant')),...
        num2str(round(imgStats.av(iav).lateWinSuppErr,2,'significant')))
    
    imgStats.av(iav).lateCycSupp = mean(mean(antiAnalysis.lateCycTC{iav}(respwin,...
        (ind)),1));
    imgStats.av(iav).lateCycSuppErr = ste(mean(antiAnalysis.lateCycTC{iav}(respwin,...
        (ind)),1),2);
    fprintf('%s Suppr. Late Cyc Resp mean/err: %s/%s\n',avName{iav},...
        num2str(round(imgStats.av(iav).lateCycSupp,2,'significant')),...
        num2str(round(imgStats.av(iav).lateCycSuppErr,2,'significant')))
    
%     x = mean(antiAnalysis.lateCycTC{visualTrials}(respwin,...
%         (cellInfo.lateCycRespCells)),1);
%     y = mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,...
%         (cellInfo.lateCycRespCells)),1);
end
[~,imgStats.lateWinTestSuppCells] = ttest(...
    mean(antiAnalysis.longTC{visualTrials}(lateWinFr,ind),1),...
    mean(antiAnalysis.longTC{auditoryTrials}(lateWinFr,ind),1));
fprintf('Late Win Test, Suppr. Cells, p=%s\n',...
    num2str(round(imgStats.lateWinTestSuppCells,2,'significant')))
[~,imgStats.lateCycTestSuppCells] = ttest(...
    mean(antiAnalysis.lateCycTC{visualTrials}(respwin,ind),1),...
    mean(antiAnalysis.lateCycTC{auditoryTrials}(respwin,ind),1));
fprintf('Late Cyc Resp Test, Suppr. Cells, p=%s\n',...
    num2str(round(imgStats.lateCycTestSuppCells,2,'significant')))

ind = antiAnalysis.lateCycAV95CITest == 1;
imgStats.modCells.nLateResp = [sum(ind' & antiAnalysis.lateCycSI' <= 0 & cellInfo.lateCycRespCells),...
    sum(ind' & antiAnalysis.lateCycSI' > 0 & cellInfo.lateCycRespCells)];
imgStats.modCells.nSupp = [sum(ind' & antiAnalysis.lateWinSI' <= 0 &  cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells)),...
    sum(ind' & antiAnalysis.lateWinSI' > 0 & cellInfo.lateSuppCells & ...
    ~(cellInfo.lateCycRespCells|cellInfo.lateRespCells|cellInfo.firstRespCells))];

fprintf('Fraction modulated late resp cells: %s\n',...
    num2str(round(sum(imgStats.modCells.nLateResp)./imgStats.nCells.lateCycResp,3,'significant')))
fprintf('Fraction modulated late supp cells: %s\n',...
    num2str(round(sum(imgStats.modCells.nSupp)./imgStats.nCells.lateSupp,3,'significant')))
fprintf('N +SI Mod. Late Resp Cells: %s/%s\n',num2str(imgStats.modCells.nLateResp(2)),...
    num2str(sum(imgStats.modCells.nLateResp)))
fprintf('N +SI Mod. Late Supp Cells: %s/%s\n',num2str(imgStats.modCells.nSupp(2)),...
    num2str(sum(imgStats.modCells.nSupp)))

save([fnout 'imgStats'],'imgStats')
%% decoding analysis figures

% setFigParams4Print('landscape')
% figure
% colormap(brewermap([],'*RdBu'));
% 
% targetRespAll = mean(targetAnalysis.tc{allTrialsInd,visualTrials}...
%     (respwin_target,:),1);
% [~,targetSortInd] = sort(targetRespAll);
% subplot 131
% % lateCycRespAll = mean(antiAnalysis.lateCycTC{allTrialsInd}(respwin,:),1);
% % [~,lateCycSortInd] = sort(lateCycRespAll);
% hm = flipud(antiAnalysis.lateCycTC{visualTrials}...
%     (tcStartFrame:cycTCEndFr,targetSortInd)');
% imagesc(hm)
% caxis(hmLim)
% colorbar
% figXAxis([],'Time (ms)',[],ttLabelFr_cyc,ttLabel_cyc)
% figYAxis([],'Cell # (Target Resp Sorted)',[])
% figAxForm
% title('Late Dist. Resp, All Cells')
% 
% subplot 132
% % targetRespAll = mean(targetAnalysis.tc{allTrialsInd,visualTrials}...
% %     (respwin_target,:),1);
% hm = flipud(targetAnalysis.tc{1,visualTrials}...
%     ((tcStartFrame-1):(cycTCEndFr-1),targetSortInd)');
% imagesc(hm)
% caxis(hmLim)
% colorbar
% figXAxis([],'Time (ms)',[],ttLabelFr_cyc,ttLabel_cyc)
% figYAxis([],'Cell # (Target Resp Sorted)',[])
% figAxForm
% title('Hard Target Resp, All Cells')
% 
% subplot 133
% hm = flipud(targetAnalysis.tc{2,visualTrials}...
%     ((tcStartFrame-1):(cycTCEndFr-1),targetSortInd)');
% imagesc(hm)
% caxis(hmLim)
% colorbar
% figXAxis([],'Time (ms)',[],ttLabelFr_cyc-1,ttLabel_cyc)
% figYAxis([],'Cell # (Target Resp Sorted)',[])
% figAxForm
% title('Easy Target Resp, All Cells')
% 
% print([fnout 'heatmapVisTrialsAllCells_dist&target'],'-dpdf','-fillpage')
% 
% figure
% colormap gray
% suptitle('Responsive Cells IDs, sorted by target response')
% 
% subplot 131
% ind = cellInfo.lateCycRespCells;
% indBW = ind; indBW(ind) = 0; indBW(~ind) = 1;
% imagesc(flipud(indBW(targetSortInd)))
% title(sprintf('Late Cyc. Resp. Cells n=%s',num2str(sum(ind))))
% figYAxis([],'Cell # (Target Resp Sorted)',[])
% figAxForm
% 
% subplot 132
% ind = cellInfo.targetRespCells;
% indBW = ind; indBW(ind) = 0; indBW(~ind) = 1;
% imagesc(flipud(indBW(targetSortInd)))
% title(sprintf('Target Resp. Cells n=%s',num2str(sum(ind))))
% figYAxis([],'Cell # (Target Resp Sorted)',[])
% figAxForm
% 
% print([fnout 'cellIDforHeatmap_dist&target'],'-dpdf','-fillpage')

if strcmp(ds,'FSAV_attentionV1')
   
    dcModel = struct;
    dcModel(1).cellTypeName = {'Tar.';'Dist.';'N.P.'};
    dcModel(1).name = 'Target';
    dcModel(2).name = 'Detect';
    cellTypeIDAllCells = sum(cat(2,(cellInfo.lateCycRespCells & ~cellInfo.targetRespCells).*1,...
        (cellInfo.lateCycRespCells & cellInfo.targetRespCells).*2,...
        (~cellInfo.lateCycRespCells & cellInfo.targetRespCells).*3),2);
    dcModel(1).cellTypeID = cellTypeIDAllCells(cellInfo.dcModelCells);
    for imod = 1:2
        dcModel(imod).av(1).name = 'Visual';
        dcModel(imod).av(2).name = 'Auditory';
        for iav = 1:2
            dcModel(2).av(iav).pctYes = nan(3,nexp);
            dcModel(imod).av(iav).dv = nan(1,nexp);
            dcModel(imod).av(iav).correlation = [];
            dcModel(imod).av(iav).weight = [];
            dcModel(imod).av(iav).pctCorrect = nan(4,nexp);
            dcModel(imod).av(iav).testPerformance = nan(4,nexp);
            dcModel(imod).av(iav).pctCorrect_movRespWin = nan(nMovWin,nexp);
            dcModel(imod).av(iav).pctCorrect_otherAV = nan(1,nexp);
        end
    end
    for iav = 1:2
        for iexp = 1:nexp
            d = decodeAnalysis(iexp).av(iav);
            dcModel(1).av(iav).dv(iexp) = d.dvTarget;
            dcModel(2).av(iav).dv(iexp) = d.dvDetect;
%             dcModel(2).av(iav).pctYes(:,iexp) = d.pctYes;

            dcModel(1).av(iav).correlation = cat(2, ...
                dcModel(1).av(iav).correlation,d.correlationTarget);
            dcModel(1).av(iav).weight = cat(1, ...
                dcModel(1).av(iav).weight,d.weightTarget);
            dcModel(1).av(iav).pctCorrect(:,iexp) = cat(2,...
                d.pctCorrectXStimTarget_holdout,d.pctCorrectAllTarget_holdout);
            dcModel(1).av(iav).testPerformance(:,iexp) = cat(2,...
                d.pctCorrectXStimTarget_train,d.pctCorrectAllTarget_train) - ...
                cat(2,...
                d.pctCorrectXStimTarget_holdout,d.pctCorrectAllTarget_holdout);
            dcModel(1).av(iav).pctCorrect_movRespWin(:,iexp) = d.pctCorrectTargetMovRespWin_holdout;
            dcModel(1).av(iav).pctCorrect_otherAV(iexp) = d.pctCorrectTarget_otherAV;

            dcModel(2).av(iav).correlation = cat(2, ...
                dcModel(2).av(iav).correlation,d.correlationDetect);
            dcModel(2).av(iav).weight = cat(1, ...
                dcModel(2).av(iav).weight,d.weightDetect);
            dcModel(2).av(iav).pctCorrect(:,iexp) = cat(2,...
                d.pctCorrectXStimDetect_holdout,d.pctCorrectAllDetect_holdout);
            dcModel(2).av(iav).testPerformance(:,iexp) = cat(2,...
                d.pctCorrectXStimDetect_train,d.pctCorrectAllDetect_train) - ...
                cat(2,...
                d.pctCorrectXStimDetect_holdout,d.pctCorrectAllDetect_holdout);
            dcModel(2).av(iav).pctCorrect_movRespWin(:,iexp) = d.pctCorrectDetectMovRespWin_holdout;
            dcModel(2).av(iav).pctCorrect_otherAV(iexp) = d.pctCorrectDetect_otherAV;
        end
    end
    stimLabel = {'0';'HT';'ET';'All'};
    mdlXStimAlpha = 0.05./3;

% response matrix (Figure 4 Schematic)

iexp=10;
for iav = 1:2
    dc = decodeDataExpt(iexp).av(iav);
    ind = decodeAnalysis(iexp).cellInd;
%     trialInd = dcTrialsUsed{iexp,iav};
%     trOut = dc.outcome(trialInd);
    trOut = dc.outcome;

%     r = dc.tc(:,ind,trialInd);
    r = dc.tc(:,ind,:);
    [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut);
    wDetect = decodeAnalysis(iexp).av(iav).weightDetect;
    wTarget = decodeAnalysis(iexp).av(iav).weightTarget;
    [~,detectWeightSort] = sort(wDetect);
    [~,targetWeightSort] = sort(wTarget);

    rDetect = {mean(r(:,:,detectTrInd==0),3);mean(r(:,:,detectTrInd==1),3)};
    rTarget = {mean(r(:,:,targetTrInd==0),3);mean(r(:,:,targetTrInd==1),3)};
    rDetectErr = {ste(r(:,:,detectTrInd==0),3);ste(r(:,:,detectTrInd==1),3)};
    rTargetErr = {ste(r(:,:,targetTrInd==0),3);ste(r(:,:,targetTrInd==1),3)};

    dfig_n = figure;
    suptitle(sprintf('%s, %s, No Resp, Detect Weight Sorted',...
        decodeDataExpt(iexp).exptName,avName{iav}))
    dfig_y = figure;
    suptitle(sprintf('%s, %s, Yes Resp, Detect Weight Sorted',...
        decodeDataExpt(iexp).exptName,avName{iav}))
    tfig_d = figure;
    suptitle(sprintf('%s, %s, Dist. Resp, Target Weight Sorted',...
        decodeDataExpt(iexp).exptName,avName{iav}))
    tfig_t = figure;
    suptitle(sprintf('%s, %s, Tar. Resp, Target Weight Sorted',...
        decodeDataExpt(iexp).exptName,avName{iav}))
    for icell = 1:sum(ind)
        cInd = detectWeightSort(icell);
        figure(dfig_n);
        subplot(5,2,icell)
%         plot(tt_cycTC,rDetect{1}(tcStartFrame:nFr_cyc,cInd),'k-')
        shadedErrorBar_chooseColor(tt_cycTC,rDetect{1}(tcStartFrame:nFr_cyc,cInd),...
            rDetectErr{1}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
        figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
        figYAxis([],'dF/F',cellRespTCLim)  
        hline(0,'r:')
        vline(respWinTT,'k--')
        figAxForm
        title(sprintf('N Resp., Cell #%s, W=%s',num2str(cInd),...
            num2str(round(wDetect(cInd),2,'significant'))))
        figure(dfig_y);
        subplot(5,2,icell)
%         plot(tt_cycTC,rDetect{2}(tcStartFrame:nFr_cyc,cInd),'k-')
        shadedErrorBar_chooseColor(tt_cycTC,rDetect{2}(tcStartFrame:nFr_cyc,cInd),...
            rDetectErr{2}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
        figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
        figYAxis([],'dF/F',cellRespTCLim)  
        hline(0,'r:')
        vline(respWinTT,'k--')
        figAxForm
        title(sprintf('Y Resp., Cell #%s, W=%s',num2str(cInd),...
            num2str(round(wDetect(cInd),2,'significant'))))
    end
    figure(dfig_n);
    print(...
        [fnout 'detectModel_' avName{iav}(1:3) '_noResp_' decodeDataExpt(iexp).exptName],...
        '-dpdf','-fillpage')
    figure(dfig_y);
    print(...
        [fnout 'detectModel_' avName{iav}(1:3) '_yesResp_' decodeDataExpt(iexp).exptName],...
        '-dpdf','-fillpage')
    for icell = 1:sum(ind)
        figure(tfig_d);
        subplot(5,2,icell)
        cInd = targetWeightSort(icell);
%         plot(tt_cycTC,rTarget{1}(tcStartFrame:nFr_cyc,cInd),'k-')
        shadedErrorBar_chooseColor(tt_cycTC,rTarget{1}(tcStartFrame:nFr_cyc,cInd),...
            rTargetErr{1}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
        figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
        figYAxis([],'dF/F',cellRespTCLim)  
        hline(0,'r:')
        vline(respWinTT,'k--')
        figAxForm
        title(sprintf('D Resp., Cell #%s, W=%s',num2str(cInd),...
            num2str(round(wTarget(cInd),2,'significant'))))
        figure(tfig_t);
        subplot(5,2,icell)
%         plot(tt_cycTC,rTarget{2}(tcStartFrame:nFr_cyc,cInd),'k-')
        shadedErrorBar_chooseColor(tt_cycTC,rTarget{2}(tcStartFrame:nFr_cyc,cInd),...
            rTargetErr{2}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
        figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
        figYAxis([],'dF/F',cellRespTCLim)  
        hline(0,'r:')
        vline(respWinTT,'k--')
        figAxForm
        title(sprintf('T Resp., Cell #%s, W=%s',num2str(cInd),...
            num2str(round(wTarget(cInd),2,'significant'))))
    end
    figure(tfig_d);
    print(...
        [fnout 'targetModel_' avName{iav}(1:3) '_distResp_' decodeDataExpt(iexp).exptName],...
        '-dpdf','-fillpage')
    figure(tfig_t);
    print(...
        [fnout 'targetModel_' avName{iav}(1:3) '_tarResp_' decodeDataExpt(iexp).exptName],...
        '-dpdf','-fillpage')
end
    
% within modality decoding in V1  (Figure 4)
    figure;
    for imod = 1:2
        if imod == 1
            iplot = 0;
        else
            iplot = 2;
        end
        for iav = 1:2
            y = dcModel(imod).av(iav).pctCorrect;
            yerr = ste(dcModel(imod).av(iav).pctCorrect,2);
            stimInd = sum(~isnan(y),2) > minTrN;
            y = y(stimInd,:);
            yerr = yerr(stimInd,:);
            x = 1:sum(stimInd);
            [mdlTest,mdlP] = ttest(y,0.5,'dim',2,'alpha',mdlXStimAlpha);
            subplot(2,2,iplot+iav)
            for iexp = 1:nexp
    %             ind = ~isnan(y(:,iexp));
                hold on
                plot(1:length(x),y(:,iexp),'k.-');
            end
            for itest = 1:length(x)
                if mdlTest(itest) == 1 
                    text(x(itest),1,sprintf('*%s',...
                        num2str(round(mdlP(itest),2,'significant'))))
                else
                    text(x(itest),1,sprintf('%s',...
                        num2str(round(mdlP(itest),2,'significant'))))
                end
            end
            errorbar(1:length(x),nanmean(y,2),yerr,'.-','MarkerSize',15)
            figXAxis([],'Stim',[0 length(x)+1],x,stimLabel(stimInd))
            hline(0.5,'k:')
            figYAxis([],'Fraction Correct',[0 1])
            figAxForm
            title(sprintf('%s Model, Train %s, Test %s',dcModel(imod).name,...
                dcModel(imod).av(iav).name,dcModel(imod).av(iav).name))
        end
    end

    print([fnout 'decodeModelPctCorrectXStim'],'-dpdf','-fillpage')

    figure;
    suptitle('Performace of trained model on hold-out trials')
    for imod = 1:2
        if imod == 1
            iplot = 0;
        else
            iplot = 2;
        end
        for iav = 1:2
            y = dcModel(imod).av(iav).testPerformance;
            yerr = ste(dcModel(imod).av(iav).testPerformance,2);
            stimInd = sum(~isnan(y),2) > minTrN;
            y = y(stimInd,:);
            yerr = yerr(stimInd,:);
            x = 1:sum(stimInd);
            [mdlTest,mdlP] = ttest(y,0,'dim',2,'alpha',mdlXStimAlpha);
            subplot(2,2,iplot+iav)
            for iexp = 1:nexp
                ind = ~isnan(y(:,iexp));
                hold on
                plot(x(ind),y(ind,iexp),'k.-');
            end
            for itest = 1:length(x)
                if mdlTest(itest) == 1
                    text(x(itest),0.5,sprintf('*%s',...
                        num2str(round(mdlP(itest),2,'significant'))))
                else
                    text(x(itest),0.5,sprintf('%s',...
                        num2str(round(mdlP(itest),2,'significant'))))
                end
            end
            errorbar(x,nanmean(y,2),yerr,'.-','MarkerSize',15)
            figXAxis([],'Stim',[0 length(x)+1],1:length(x),stimLabel(stimInd))
            hline(0,'k:')
            figYAxis([],'Train - HO Frac. Correct',[-0.5 0.5])
            figAxForm
            title(sprintf('%s Model, Train %s, Test %s',dcModel(imod).name,...
                dcModel(imod).av(iav).name,dcModel(imod).av(iav).name))
        end
    end

    print([fnout 'decodeModelPerformanceXStim'],'-dpdf','-fillpage')

% timing of decoding (Suppl. Fig. 8)
    figure
    suptitle('Percent Correct of responses analyzed at different time windows around stimulus; red is used analysis window; blue is min RT')
    winInd = [1:3,5:nwins];
    x = movWinLabelMs(winInd);
    pctCorrThresh = 0.55;
    for imod = 1:2
        if imod == 1
            iplot = 0;
        else
            iplot = 2;
        end
        for iav = 1:2
            y = dcModel(imod).av(iav).pctCorrect_movRespWin(winInd,:);
            yerr = ste(dcModel(imod).av(iav).pctCorrect_movRespWin(winInd,:),2);
            targetTimes = x(x>0);
            firstTimes = nan(1,length(yerr));
            for i = 1:length(yerr)
                ind = find(y(x>0,i)>pctCorrThresh,1);
                firstTimes(i) = targetTimes(ind);
            end
            subplot(2,2,iplot+iav)
            for iexp = 1:nexp
                ind = ~isnan(y(:,iexp));
                hold on
                plot(x(ind),y(ind,iexp),'k.-');
            end
            errorbar(x,nanmean(y,2),yerr,'.-','MarkerSize',15)
            errorbar(mean(firstTimes),0.9,[],[],ste(firstTimes,2),ste(firstTimes,2),'.')
            text(0,0.95,[num2str(round(mean(firstTimes),2,'significant'))...
                '+/-' num2str(round(ste(firstTimes,2),2,'significant'))])
            figXAxis([],'Resp. Win. Center (Relative to Stim)',[x(1) x(end)],x(1:2:length(x)),round(x(1:2:length(x))))
            hline(0.5,'k:')
            figYAxis([],'Frac. Correct',[0 1])
            figAxForm
            title(sprintf('%s Model, Train %s, Test %s',dcModel(imod).name,...
                dcModel(imod).av(iav).name,dcModel(imod).av(iav).name))
            vline(respWinTT,'r--')
            vline(baseWinTT,'k--')
            vline(0,'k:')
            if iav == 1
                vline(visRTwindow(1),'b--')
            else
                vline(audRTwindow(1),'b--')
            end
        end
    end

    print([fnout 'decodeModelPerformanceXRespWindow'],'-dpdf','-fillpage')

% within modality weights    
for iav = 1:2
        figure
        suptitle(sprintf('%s,Task correlation and weights of cells used in model',...
            dcModel(1).av(iav).name))

        weightInd = dcModel(imod).av(iav).weight > 0;
        binnedTargetWeight = nan(1,2);
        binnedTargetWeightErr = nan(1,2);
        binnedDetectWeight = nan(1,2);
        binnedDetectWeightErr = nan(1,2);
        binnedDetectWeightTest = nan(1,2);
        binnedTargetWeight(1) = mean(dcModel(1).av(iav).weight(~weightInd));
        binnedTargetWeight(2) = mean(dcModel(1).av(iav).weight(weightInd));
        binnedTargetWeightErr(1) = ste(dcModel(1).av(iav).weight(~weightInd),1);
        binnedTargetWeightErr(2) = ste(dcModel(1).av(iav).weight(weightInd),1);
        binnedDetectWeight(1) = mean(dcModel(2).av(iav).weight(~weightInd));
        binnedDetectWeight(2) = mean(dcModel(2).av(iav).weight(weightInd));
        binnedDetectWeightErr(1) = ste(dcModel(2).av(iav).weight(~weightInd),1);
        binnedDetectWeightErr(2) = ste(dcModel(2).av(iav).weight(weightInd),1);
        [~,binnedDetectWeightTest(1)] = ttest(dcModel(2).av(iav).weight(~weightInd));
        [~,binnedDetectWeightTest(2)] = ttest(dcModel(2).av(iav).weight(weightInd));

        for imod = 1:2
            x = dcModel(imod).av(iav).correlation;
            y = dcModel(imod).av(iav).weight;

            subplot(2,2,imod)
            plot(x,y,'.')
            figXAxis([],'Correlation',[-1 1])
            figYAxis([],'Weight',weightLim)
            figAxForm
            vline(0,'k:')
            hline(0,'k:')
            title(sprintf('%s Model',dcModel(imod).name))
        end

        subplot(2,2,3)
        x = dcModel(1).av(iav).weight;
        y = dcModel(2).av(iav).weight;
        [corrmat,pmat] = corrcoef(x,y);
        plot(x,y,'.')
        figXAxis([],'Target Weight',[-2 3])
        figYAxis([],'Detect Weight',[-2 3])
        figAxForm
        hline(0,'k:')
        vline(0,'k:')
        title(sprintf('Corr=%s, p=%s',num2str(round(corrmat(1,2),2,'significant')),...
            num2str(round(pmat(1,2),2,'significant'))))

        subplot(2,2,4)
        hold on
        for i = 1:2
            if i == 1
                ind = ~weightInd;
            elseif i == 2
                ind = weightInd;
            end
            y = dcModel(2).av(iav).weight(ind);
            x = binnedTargetWeight(i).*(ones(sum(ind),1));
            plot(x,y,'.')
            text(binnedTargetWeight(i),weightLim(end),num2str(round(...
                binnedDetectWeightTest(i),2,'significant')))
        end
        errorbar(binnedTargetWeight,binnedDetectWeight,binnedDetectWeightErr,'.-')
        figXAxis([],'Target Weight',binnedWeightLim)
        figYAxis([],'Detect Weight',weightLim)
        figAxForm
        hline(0,'k:')
        vline(0,'k:')

        dcModel(1).av(iav).binnedWeight = binnedTargetWeight;
        dcModel(1).av(iav).binnedWeightErr = binnedTargetWeightErr;
        dcModel(2).av(iav).binnedWeight = binnedDetectWeight;
        dcModel(2).av(iav).binnedWeightErr = binnedDetectWeightErr;
        dcModel(2).av(iav).binnedWeightTest = binnedDetectWeightTest;

        print([fnout 'decodeModelWeights_' dcModel(1).av(iav).name],'-dpdf','-fillpage')
end

% orthogonal coding (Figure 5)
     
    WxAttnFig = figure;
    suptitle('Target Model Attention - Late Resp. Cells')
    WFig = figure;
    suptitle('Target Model Weight x vis AuROC Group - Late Resp. Cells')
    weightDataFor2way = cell(1,2);
    weightIndFor2way = cell(1,2);
    for iav = 1:2
        figure(WxAttnFig)
        subplot(2,2,iav)
    %     ind = cellInfo.lateCycRespCells;
        ind1 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.lateStimAuROCTest...
            & cellInfo.lateStimAuROC > 0.5;
        ind2 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.lateStimAuROCTest...
            & cellInfo.lateStimAuROC < 0.5;
        ind3 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & ~cellInfo.lateStimAuROCTest;
        n = {[dcModel(1).cellTypeName{1} ';n=' num2str(sum(ind1))];...
            [dcModel(1).cellTypeName{2} ';n=' num2str(sum(ind2))];...
            [dcModel(1).cellTypeName{3} ';n=' num2str(sum(ind3))]};

        x = antiAnalysis.lateCycSI(ind1);
        y = cellInfo.targetWeight{iav}(ind1);    
        plot(x,y,'.')
        hold on
        x = antiAnalysis.lateCycSI(ind2);
        y = cellInfo.targetWeight{iav}(ind2);    
        plot(x,y,'.')
        x = antiAnalysis.lateCycSI(ind3);
        y = cellInfo.targetWeight{iav}(ind3);    
        plot(x,y,'.')
        [corrmat,pmat] = corrcoef(antiAnalysis.lateCycSI(cellInfo.lateCycRespCells),...
            cellInfo.targetWeight{iav}(cellInfo.lateCycRespCells));
        r = corrmat(1,2);
        p = pmat(1,2);

        figXAxis([],'Selectivity',siLim)
        figYAxis([],'Weight',weightLim)
        figAxForm
        vline(0,'k:')
        hline(0,'k:')
        title(sprintf('Train %s, R=%s, p=%s',dcModel(1).av(iav).name,...
            num2str(round(r,2,'significant')),num2str(round(p,2,'significant'))))

        subplot(2,2,iav+2)

        x = mean(antiAnalysis.lateCycSI(ind1));
        xerr = ste(antiAnalysis.lateCycSI(ind1),2);
        y = nanmean(cellInfo.targetWeight{iav}(ind1));
        yerr = ste(cellInfo.targetWeight{iav}(ind1),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')
        hold on
        x = mean(antiAnalysis.lateCycSI(ind2));
        xerr = ste(antiAnalysis.lateCycSI(ind2),2);
        y = nanmean(cellInfo.targetWeight{iav}(ind2));
        yerr = ste(cellInfo.targetWeight{iav}(ind2),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')
        x = mean(antiAnalysis.lateCycSI(ind3));
        xerr = ste(antiAnalysis.lateCycSI(ind3),2);
        y = nanmean(cellInfo.targetWeight{iav}(ind3));
        yerr = ste(cellInfo.targetWeight{iav}(ind3),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')

        figXAxis([],'Selectivity',siLimSum)
        figYAxis([],'Weight',weightLimSum)
        figAxForm
        vline(0,'k:')
        hline(0,'k:')
        title(sprintf('Train %s',dcModel(1).av(iav).name))
        L = legend(n,'location','southwest');
        title(L,'Neuron Pref. in Vis Trials')
        
        figure(WFig)
        pt = [];
        subplot(1,2,iav)
        y = nanmean(cellInfo.targetWeight{iav}(ind1));
        yerr = ste(cellInfo.targetWeight{iav}(ind1),2);    
        errorbar(1,y,yerr,yerr,'.')
        [~,pt(1)] = ttest(cellInfo.targetWeight{iav}(ind1));
        y = nanmean(cellInfo.targetWeight{iav}(ind2));
        yerr = ste(cellInfo.targetWeight{iav}(ind2),2);   
        hold on
        errorbar(2,y,yerr,yerr,'.')
        [~,pt(2)] = ttest(cellInfo.targetWeight{iav}(ind2));
        y = nanmean(cellInfo.targetWeight{iav}(ind3));
        yerr = ste(cellInfo.targetWeight{iav}(ind3),2);    
        errorbar(3,y,yerr,yerr,'.')
        [~,pt(3)] = ttest(cellInfo.targetWeight{iav}(ind3));
        x = cat(2,ones(1,sum(ind1)),ones(1,sum(ind2)).*2,ones(1,sum(ind3)).*3);
        y = cat(2,cellInfo.targetWeight{iav}(ind1),...
            cellInfo.targetWeight{iav}(ind2),cellInfo.targetWeight{iav}(ind3));
        [p,~,stats] = anova1(y,x,'off');
        posthoc = multcompare(stats,[],'off');
        fprintf('%s Target Weight x Vis auROC Model\n',avName{iav})
        disp(posthoc(:,[1,2,end]))
        fprintf('%s Target Weight Vs. Zero x Vis auROC\n',avName{iav})
        disp(round(pt,2,'significant'))
        figXAxis([],'Vis. auROC Group',[0 4],1:3,{'T','D','NP'})
        figYAxis([],'Weight',weightLimSum)
        figAxForm
        hline(0,'k:')
        title(sprintf('Train %s, p=%s',dcModel(1).av(iav).name,...
            num2str(round(p,2,'significant'))))
        L = legend(n,'location','southwest');
        title(L,'Neuron Pref. in Vis Trials, p=%s')
        
        weightDataFor2way{iav} = y;
        weightIndFor2way{iav} = x;
    end
%     [p,~,stats] = anovan(cell2mat(weightDataFor2way),...
%         {cell2mat(weightIndFor2way) cat(2,ones(1,length(weightDataFor2way{1})),...
%         ones(1,length(weightDataFor2way{1}))*2)},'model','interaction');
%     multcompare(stats)
    
    figure(WFig)
    print([fnout 'targetModelWeightxGroup'],'-dpdf')
    figure(WxAttnFig)
    print([fnout 'targetModelWeightXAttn'],'-dpdf','-fillpage')

    WxAttnFig = figure;
    suptitle('Detect Model Attention - Late Resp. Cells')
    WFig = figure;
    suptitle('Detect Model Weight x vis AuROC Group - Late Resp. Cells')
    for iav = 1:2
        figure(WxAttnFig)
        subplot(2,2,iav)
    %     ind = cellInfo.lateCycRespCells;
        ind1 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.lateStimAuROCTest...
            & cellInfo.lateStimAuROC > 0.5;
        ind2 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.lateStimAuROCTest...
            & cellInfo.lateStimAuROC < 0.5;
        ind3 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & ~cellInfo.lateStimAuROCTest;
        n = {[dcModel(1).cellTypeName{1} ';n=' num2str(sum(ind1))];...
            [dcModel(1).cellTypeName{2} ';n=' num2str(sum(ind2))];...
            [dcModel(1).cellTypeName{3} ';n=' num2str(sum(ind3))]};

        x = antiAnalysis.lateCycSI(ind1);
        y = cellInfo.detectWeight{iav}(ind1);    
        plot(x,y,'.')
        hold on
        x = antiAnalysis.lateCycSI(ind2);
        y = cellInfo.detectWeight{iav}(ind2);    
        plot(x,y,'.')
        x = antiAnalysis.lateCycSI(ind3);
        y = cellInfo.detectWeight{iav}(ind3);    
        plot(x,y,'.')
        [corrmat,pmat] = corrcoef(antiAnalysis.lateCycSI(cellInfo.lateCycRespCells),...
            cellInfo.detectWeight{iav}(cellInfo.lateCycRespCells));
        r = corrmat(1,2);
        p = pmat(1,2);

        figXAxis([],'Selectivity',siLim)
        figYAxis([],'Weight',weightLim)
        figAxForm
        vline(0,'k:')
        hline(0,'k:')
        title(sprintf('Train %s, R=%s, p=%s',dcModel(2).av(iav).name,...
            num2str(round(r,2,'significant')),num2str(round(p,2,'significant'))))

        subplot(2,2,iav+2)

        x = mean(antiAnalysis.lateCycSI(ind1));
        xerr = ste(antiAnalysis.lateCycSI(ind1),2);
        y = nanmean(cellInfo.detectWeight{iav}(ind1));
        yerr = ste(cellInfo.detectWeight{iav}(ind1),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')
        hold on
        x = mean(antiAnalysis.lateCycSI(ind2));
        xerr = ste(antiAnalysis.lateCycSI(ind2),2);
        y = nanmean(cellInfo.detectWeight{iav}(ind2));
        yerr = ste(cellInfo.detectWeight{iav}(ind2),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')
        x = mean(antiAnalysis.lateCycSI(ind3));
        xerr = ste(antiAnalysis.lateCycSI(ind3),2);
        y = nanmean(cellInfo.detectWeight{iav}(ind3));
        yerr = ste(cellInfo.detectWeight{iav}(ind3),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')

        figXAxis([],'Selectivity',siLimSum)
        figYAxis([],'Weight',weightLimSum)
        figAxForm
        vline(0,'k:')
        hline(0,'k:')
        title(sprintf('Train %s',dcModel(1).av(iav).name))
        L = legend(n,'location','southwest');
        title(L,'Neuron Pref. in Vis Trials')
        
        figure(WFig)
        subplot(1,2,iav)
        y = nanmean(cellInfo.detectWeight{iav}(ind1));
        yerr = ste(cellInfo.detectWeight{iav}(ind1),2);    
        errorbar(1,y,yerr,yerr,'.')
        [~,pt(1)] = ttest(cellInfo.targetWeight{iav}(ind1));
        y = nanmean(cellInfo.detectWeight{iav}(ind2));
        yerr = ste(cellInfo.detectWeight{iav}(ind2),2);   
        hold on
        errorbar(2,y,yerr,yerr,'.')
        [~,pt(2)] = ttest(cellInfo.targetWeight{iav}(ind2));
        y = nanmean(cellInfo.detectWeight{iav}(ind3));
        yerr = ste(cellInfo.detectWeight{iav}(ind3),2);    
        errorbar(3,y,yerr,yerr,'.')
        [~,pt(3)] = ttest(cellInfo.targetWeight{iav}(ind3));
        x = cat(2,ones(1,sum(ind1)),ones(1,sum(ind2)).*2,ones(1,sum(ind3)).*3);
        y = cat(2,cellInfo.detectWeight{iav}(ind1),...
            cellInfo.detectWeight{iav}(ind2),cellInfo.detectWeight{iav}(ind3));
        [p,~,stats] = anova1(y,x,'off');
        posthoc = multcompare(stats,[],'off');
        fprintf('%s Detect Weight x Vis auROC Model\n',avName{iav})
        disp(posthoc(:,[1,2,end]))
        fprintf('%s Target Weight Vs. Zero x Vis auROC\n',avName{iav})
        disp(round(pt,2,'significant'))
        figXAxis([],'Vis. auROC Group',[0 4],1:3,{'T','D','NP'})
        figYAxis([],'Weight',weightLimSum)
        figAxForm
        hline(0,'k:')
        title(sprintf('Train %s,p=%s',dcModel(1).av(iav).name,...
            num2str(round(p,2,'significant'))))
        L = legend(n,'location','southwest');
        title(L,'Neuron Pref. in Vis Trials')
    end
    
    figure(WFig)
    print([fnout 'detectModelWeightxGroup'],'-dpdf')
    figure(WxAttnFig)
    print([fnout 'detectModelWeightXAttn'],'-dpdf','-fillpage')

    WxAttnFig = figure;
    suptitle('Target Model Attention - Late Resp. Cells')
    WFig = figure;
    suptitle('Target Model Weight x aud AuROC Group - Late Resp. Cells')
    for iav = 1:2
        figure(WxAttnFig)
        subplot(2,2,iav)
    %     ind = (cellInfo.lateCycRespCells|cellInfo.targetRespCells);
        ind1 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.audLateStimAuROCTest...
            & cellInfo.audLateStimAuROC > 0.5;
        ind2 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.audLateStimAuROCTest...
            & cellInfo.audLateStimAuROC < 0.5;
        ind3 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & ~cellInfo.audLateStimAuROCTest;
        n = {[dcModel(1).cellTypeName{1} ';n=' num2str(sum(ind1))];...
            [dcModel(1).cellTypeName{2} ';n=' num2str(sum(ind2))];...
            [dcModel(1).cellTypeName{3} ';n=' num2str(sum(ind3))]};

        x = antiAnalysis.lateCycSI(ind1);
        y = cellInfo.targetWeight{iav}(ind1);    
        plot(x,y,'.')
        hold on
        x = antiAnalysis.lateCycSI(ind2);
        y = cellInfo.targetWeight{iav}(ind2);    
        plot(x,y,'.')
        x = antiAnalysis.lateCycSI(ind3);
        y = cellInfo.targetWeight{iav}(ind3);    
        plot(x,y,'.')

        figXAxis([],'Selectivity',siLim)
        figYAxis([],'Weight',weightLim)
        figAxForm
        vline(0,'k:')
        hline(0,'k:')
        title(sprintf('Train %s',dcModel(1).av(iav).name))

        subplot(2,2,iav+2)

        x = mean(antiAnalysis.lateCycSI(ind1));
        xerr = ste(antiAnalysis.lateCycSI(ind1),2);
        y = nanmean(cellInfo.targetWeight{iav}(ind1));
        yerr = ste(cellInfo.targetWeight{iav}(ind1),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')
        hold on
        x = mean(antiAnalysis.lateCycSI(ind2));
        xerr = ste(antiAnalysis.lateCycSI(ind2),2);
        y = nanmean(cellInfo.targetWeight{iav}(ind2));
        yerr = ste(cellInfo.targetWeight{iav}(ind2),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')
        x = mean(antiAnalysis.lateCycSI(ind3));
        xerr = ste(antiAnalysis.lateCycSI(ind3),2);
        y = nanmean(cellInfo.targetWeight{iav}(ind3));
        yerr = ste(cellInfo.targetWeight{iav}(ind3),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')

        figure(WFig)
        subplot(1,2,iav)
        y = nanmean(cellInfo.targetWeight{iav}(ind1));
        yerr = ste(cellInfo.targetWeight{iav}(ind1),2);    
        errorbar(1,y,yerr,yerr,'.')
        [~,pt(1)] = ttest(cellInfo.targetWeight{iav}(ind1));
        y = nanmean(cellInfo.targetWeight{iav}(ind2));
        yerr = ste(cellInfo.targetWeight{iav}(ind2),2);   
        hold on
        errorbar(2,y,yerr,yerr,'.')
        [~,pt(2)] = ttest(cellInfo.targetWeight{iav}(ind2));
        y = nanmean(cellInfo.targetWeight{iav}(ind3));
        yerr = ste(cellInfo.targetWeight{iav}(ind3),2);    
        errorbar(3,y,yerr,yerr,'.')
        [~,pt(3)] = ttest(cellInfo.targetWeight{iav}(ind3));
        x = cat(2,ones(1,sum(ind1)),ones(1,sum(ind2)).*2,ones(1,sum(ind3)).*3);
        y = cat(2,cellInfo.targetWeight{iav}(ind1),...
            cellInfo.targetWeight{iav}(ind2),cellInfo.targetWeight{iav}(ind3));
        [p,~,stats] = anova1(y,x,'off');
        posthoc = multcompare(stats,[],'off');
        fprintf('%s Target Weight x Aud auROC Model',avName{iav})
        disp(posthoc(:,[1,2,end]))
        fprintf('%s Target Weight Vs. Zero x Aud auROC\n',avName{iav})
        disp(round(pt,2,'significant'))
        figXAxis([],'auROC Group',[0 4],1:3,{'T','D','NP'})
        figYAxis([],'Weight',weightLimSum)
        figAxForm
        hline(0,'k:')
        title(sprintf('Train %s,p=%s',dcModel(1).av(iav).name,...
            num2str(round(p,2,'significant'))))
        L = legend(n,'location','southwest');
        title(L,'Neuron Pref. in Aud Trials')
    end

    figure(WFig)
    print([fnout 'targetModelAudSortWeightXGroup'],'-dpdf')
    figure(WxAttnFig)
    print([fnout 'targetModelAudSortWeightXAttn'],'-dpdf','-fillpage')

    WxAttnFig = figure;
    suptitle('Detect Model Attention - Late Resp. Cells')
    WFig = figure;
    suptitle('Detect Model Weight x aud AuROC Group - Late Resp. Cells')
    for iav = 1:2
        figure(WxAttnFig)
        subplot(2,2,iav)
    %     ind = (cellInfo.lateCycRespCells|cellInfo.targetRespCells);
        ind1 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.audLateStimAuROCTest...
            & cellInfo.audLateStimAuROC > 0.5;
        ind2 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & cellInfo.audLateStimAuROCTest...
            & cellInfo.audLateStimAuROC < 0.5;
        ind3 = (cellInfo.lateCycRespCells|cellInfo.targetRespCells) & ~cellInfo.audLateStimAuROCTest;
        n = {[dcModel(1).cellTypeName{1} ';n=' num2str(sum(ind1))];...
            [dcModel(1).cellTypeName{2} ';n=' num2str(sum(ind2))];...
            [dcModel(1).cellTypeName{3} ';n=' num2str(sum(ind3))]};

        x = antiAnalysis.lateCycSI(ind1);
        y = cellInfo.detectWeight{iav}(ind1);    
        plot(x,y,'.')
        hold on
        x = antiAnalysis.lateCycSI(ind2);
        y = cellInfo.detectWeight{iav}(ind2);    
        plot(x,y,'.')
        x = antiAnalysis.lateCycSI(ind3);
        y = cellInfo.detectWeight{iav}(ind3);    
        plot(x,y,'.')

        figXAxis([],'Selectivity',siLim)
        figYAxis([],'Weight',weightLim)
        figAxForm
        vline(0,'k:')
        hline(0,'k:')
        title(sprintf('Train %s',dcModel(1).av(iav).name))

        subplot(2,2,iav+2)

        x = mean(antiAnalysis.lateCycSI(ind1));
        xerr = ste(antiAnalysis.lateCycSI(ind1),2);
        y = nanmean(cellInfo.detectWeight{iav}(ind1));
        yerr = ste(cellInfo.detectWeight{iav}(ind1),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')
        hold on
        x = mean(antiAnalysis.lateCycSI(ind2));
        xerr = ste(antiAnalysis.lateCycSI(ind2),2);
        y = nanmean(cellInfo.detectWeight{iav}(ind2));
        yerr = ste(cellInfo.detectWeight{iav}(ind2),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')
        x = mean(antiAnalysis.lateCycSI(ind3));
        xerr = ste(antiAnalysis.lateCycSI(ind3),2);
        y = nanmean(cellInfo.detectWeight{iav}(ind3));
        yerr = ste(cellInfo.detectWeight{iav}(ind3),2);    
        errorbar(x,y,yerr,yerr,xerr,xerr,'.')

        figXAxis([],'Selectivity',siLimSum)
        figYAxis([],'Weight',weightLimSum)
        figAxForm
        vline(0,'k:')
        hline(0,'k:')
        title(sprintf('Train %s',dcModel(1).av(iav).name))
        L = legend(n,'location','southwest');
        title(L,'Neuron Pref. in Vis Trials')
        
        figure(WFig)
        subplot(1,2,iav)
        y = nanmean(cellInfo.detectWeight{iav}(ind1));
        yerr = ste(cellInfo.detectWeight{iav}(ind1),2);    
        errorbar(1,y,yerr,yerr,'.')
        [~,pt(1)] = ttest(cellInfo.targetWeight{iav}(ind1));
        y = nanmean(cellInfo.detectWeight{iav}(ind2));
        yerr = ste(cellInfo.detectWeight{iav}(ind2),2);   
        hold on
        errorbar(2,y,yerr,yerr,'.')
        [~,pt(2)] = ttest(cellInfo.targetWeight{iav}(ind2));
        y = nanmean(cellInfo.detectWeight{iav}(ind3));
        yerr = ste(cellInfo.detectWeight{iav}(ind3),2);    
        errorbar(3,y,yerr,yerr,'.')
        [~,pt(3)] = ttest(cellInfo.targetWeight{iav}(ind3));
        x = cat(2,ones(1,sum(ind1)),ones(1,sum(ind2)).*2,ones(1,sum(ind3)).*3);
        y = cat(2,cellInfo.detectWeight{iav}(ind1),...
            cellInfo.detectWeight{iav}(ind2),cellInfo.detectWeight{iav}(ind3));
        [p,~,stats] = anova1(y,x,'off');
        posthoc = multcompare(stats,[],'off');
        fprintf('%s Target Weight x Vis auROC Model',avName{iav})
        disp(posthoc(:,[1,2,end]))
        fprintf('%s Target Weight Vs. Zero x Aud auROC\n',avName{iav})
        disp(round(pt,2,'significant'))
        figXAxis([],'auROC Group',[0 4],1:3,{'T','D','NP'})
        figYAxis([],'Weight',weightLimSum)
        figAxForm
        hline(0,'k:')
        title(sprintf('Train %s,p=%s',dcModel(1).av(iav).name,...
            num2str(round(p,2,'significant'))))
        L = legend(n,'location','southwest');
        title(L,'Neuron Pref. in Aud Trials')
    end

    figure(WFig)
    print([fnout 'detectModelAudSortWeightXGroup'],'-dpdf')
    figure(WxAttnFig)
    print([fnout 'detectModelAudSortWeightXAttn'],'-dpdf','-fillpage')

    figure
    for imod = 1:2
        for iav = 1:2
            if imod == 1
                subplot(2,2,iav)
            elseif imod == 2
                subplot(2,2,iav+2)
            end
            if iav == visualTrials
                otherAV = 'Auditory';
            else
                otherAV = 'Visual';
            end
            x = 1:2;
            y = [dcModel(imod).av(iav).pctCorrect(end,:)', dcModel(imod).av(iav).pctCorrect_otherAV'];
            hold on
            plot(x,y,'k')
            errorbar(x,mean(y,1),ste(y,1),'.')
            [~,p] = ttest(y,0.5);
            [~,p2] = ttest(y(:,1),y(:,2));

            figXAxis([],'Test Trials',[0 3],x,{dcModel(imod).av(iav).name,otherAV})
            figYAxis([],'Fraction Correct',[0 1])
            figAxForm
            hline(0.5,'k:')
            title(sprintf('%s Model, Train %s',dcModel(imod).name,dcModel(imod).av(iav).name))
            for i = 1:2
                if p(i) < 0.05/7
                    text(i,1,sprintf('*%s',num2str(round(p(i),2,'significant'))))
                else
                    text(i,1,num2str(round(p(i),2,'significant')))
                end
            end
            if p2 < 0.05/4
                text(1.5,0,sprintf('*%s',num2str(round(p2,2,'significant'))))
            else
                text(1.5,0,sprintf('%s',num2str(round(p2,2,'significant'))))
            end
        end
    end

    print([fnout 'decodeModelPctCorrectOtherAVModel'],'-dpdf','-fillpage')

    % compare visual and auditory weights
    figure
    subplot 221
    ind = cellInfo.lateCycRespCells | cellInfo.targetRespCells;
    x = cellInfo.targetWeight{visualTrials}(ind)';
    y = cellInfo.targetWeight{auditoryTrials}(ind)';
    hold on
    plot(x,y,'.')
    [corrmat,pmat] = corrcoef(x,y);
    c = corrmat(1,2);
    p = pmat(1,2);
    mdl = fitlm(x,y);
    rsq = mdl.Rsquared.Ordinary;
    yfit = predict(mdl,x);
    plot(x,yfit,'-')
    [xMax,yInd] = max(x);
    plot([0 xMax],[0 y(yInd)],'-')
    figXAxis([],'Visual Weight',weightLim)
    figYAxis([],'Auditory Weight',weightLim)
    figAxForm
    title(sprintf('Target Model, Rsq = %s, corr = %s, p = %s',num2str(round(rsq,2,'significant')),...
        num2str(round(c,2,'significant')),num2str(round(p,2,'significant')))) 
    
    oriEdges = [0:15:90];
    subplot 222
    targetWeightAVangle = rad2deg(abs(atan(y./x)));
    histogram(targetWeightAVangle,oriEdges,'Normalization','probability')
    figXAxis([],'Angle of A/V Weight',[-10 100],0:15:90)
    figYAxis([],'Fraction of Cells',[0 .3])
    figAxForm
    title(sprintf('Target Model,exCell Angle=%s',num2str(round(targetWeightAVangle(yInd),2,'significant'))))
    
    
    subplot 223
    ind = cellInfo.lateCycRespCells | cellInfo.targetRespCells;
    x = cellInfo.detectWeight{visualTrials}(ind)';
    y = cellInfo.detectWeight{auditoryTrials}(ind)';
    hold on
    plot(x,y,'.')
    [corrmat,pmat] = corrcoef(x,y);
    p = pmat(1,2);
    c = corrmat(1,2);
    mdl = fitlm(x,y);
    rsq = mdl.Rsquared.Ordinary;
    yfit = predict(mdl,x);
    plot(x,yfit,'-')
    [yMax,xInd] = max(y);
    plot([0 x(xInd)],[0 yMax],'-')
    figXAxis([],'Visual Weight',weightLim)
    figYAxis([],'Auditory Weight',weightLim)
    figAxForm
    title(sprintf('Detect Model, Rsq = %s, corr = %s, p = %s',num2str(round(rsq,2,'significant')),...
        num2str(round(c,2,'significant')),num2str(round(p,2,'significant')))) 

    subplot 224
    detectWeightAVangle = rad2deg(abs(atan(y./x)));
    histogram(detectWeightAVangle,oriEdges,'Normalization','probability')
    figXAxis([],'Angle of A/V Weight',[-10 100],0:15:90)
    figYAxis([],'Fraction of Cells',[0 .3])
    figAxForm
    title(sprintf('Detect Model,exCell Angle=%s',num2str(round(detectWeightAVangle(yInd),2,'significant'))))
    
    print([fnout 'decodeModelCompareAVWeights'],'-dpdf','-fillpage')
    
    figure
    suptitle('Visual Valid and Invalid Trials Matched for Orientation')
    subplot 121
    y = [];
    for iexp = 1:nexp
        if ~isempty(decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget)
            y = cat(1,y,[decodeAnalysis(iexp).av(visualTrials).validMatchedPctCorrectTarget_holdout,...
                decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget, ...
                decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget_testAudModel]);
        end
    end
    plot(1:3,y,'k-')
    hold on
    errorbar(1:3,nanmean(y,1),ste(y,1),'.')
    [~,p]= ttest(y,0.5);
    for i = 1:3
        text(i,1,num2str(round(p(i),2,'significant')))
    end
%     [~,p] = ttest(y(:,1),y(:,2));
    [p,~,stats] = anova1(y,[],'off');
    anovaCmpTarget = multcompare(stats,'display','off');
%     text(1.5,.1,num2str(round(p,2,'significant')))
    figXAxis([],'Train-Test-Cue',[0 4],1:3,{'Vis-Vis-Val','Vis-Vis-Inv','Aud-Vis-Inv'})
    figYAxis([],'Fraction Correct',[0 1])
    figAxForm
    hline(0.5,'k:')
    title(sprintf('Target Model, p=%s',num2str(round(p,2,'significant'))))
    subplot 122
    y = [];
    for iexp = 1:nexp
        if ~isempty(decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect)
            y = cat(1,y,[decodeAnalysis(iexp).av(visualTrials).validMatchedPctCorrectDetect_holdout,...
                decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect, ...
                decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect_testAudModel]);
        end
    end
    plot(1:3,y,'k-')
    hold on
    errorbar(1:3,nanmean(y,1),ste(y,1),'.')
    [~,p]= ttest(y,0.5);
    for i = 1:3
        text(i,1,num2str(round(p(i),2,'significant')))
    end
%     [~,p] = ttest(y(:,1),y(:,2));
    [p,~,stats] = anova1(y,[],'off');
    anovaCmpDetect = multcompare(stats,'display','off');
%     text(1.5,.1,num2str(round(p,2,'significant')))
    figXAxis([],'Train-Test-Cue',[0 4],1:3,{'Vis-Vis-Val','Vis-Vis-Inv','Aud-Vis-Inv'})
    figYAxis([],'Fraction Correct',[0 1])
    figAxForm
    hline(0.5,'k:')
    title(sprintf('Detect Model, p=%s',num2str(round(p,2,'significant'))))  
    print([fnout 'decodeModelCompareInvPctCorrect'],'-dpdf')
    
% test combo model (Suppl. Fig 9)
    avName = {'Visual','Auditory'};
    setFigParams4Print('portrait')
    figure
    suptitle('p val in title is within mod vs. combo')
for iav = 1:2
    if iav == 1
        iplot = 1;
        otherAV = 2;
    else
        iplot = 3;
        otherAV = 1;
    end
    subplot(2,2,iplot)
    y = [];
    for iexp = 1:nexp
%             if iav == 1
            y = cat(1,y,[decodeAnalysis(iexp).av(iav).pctCorrectAllTarget_holdout,...
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_comboTrain]);
%             else
%                 y = cat(1,y,[decodeAnalysis(iexp).av(iav).pctCorrectTarget_audTrainComboMatch,...
%                     decodeAnalysis(iexp).av(iav).pctCorrectTarget_comboTrain]);
%             end
    end
    plot(1:2,y,'k-')
    hold on
    errorbar(1:2,mean(y,1),ste(y,1),'.')
    [~,p]= ttest(y,0.5);
    for i = 1:2
        text(i,1,num2str(round(p(i),2,'significant')))
    end
    figXAxis([],'Train Target',[0 3],1:2,{avName{iav};'Vis+Aud'})
    figYAxis([],'Fraction Correct',[0 1])
    figAxForm
    hline(0.5,'k:')
    [~,p]= ttest(y(:,1),y(:,2));
    title(sprintf('Test %s, Target Model, p=%s',avName{iav},num2str(round(p,2,'significant'))))
    subplot(2,2,iplot+1)
    y = [];
    for iexp = 1:nexp
%             if iav == 1
            y = cat(1,y,[decodeAnalysis(iexp).av(iav).pctCorrectAllDetect_holdout,...
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_comboTrain]);
%             else
%                 y = cat(1,y,[decodeAnalysis(iexp).av(iav).pctCorrectDetect_audTrainComboMatch,...
%                     decodeAnalysis(iexp).av(iav).pctCorrectTarget_comboTrain]);
%             end
    end
    plot(1:2,y,'k-')
    hold on
    errorbar(1:2,mean(y,1),ste(y,1),'.')
    [~,p]= ttest(y,0.5);
    for i = 1:2
        text(i,1,num2str(round(p(i),2,'significant')))
    end
    figXAxis([],'Train Detect',[0 3],1:2,{avName{iav};'Vis+Aud'})
    figYAxis([],'Fraction Correct',[0 1])
    figAxForm
    hline(0.5,'k:')
    [~,p]= ttest(y(:,1),y(:,2));
    title(sprintf('Test %s, Detect Model, p=%s',avName{iav},num2str(round(p,2,'significant'))))
end
    print([fnout 'decodeModelComboTrainPctCorrect'],'-dpdf','-fillpage')
    
    setFigParams4Print('landscape')
    figure;
    suptitle('Task-responsive neurons, bootstrapped weights')
    ind = cellInfo.targetRespCells | cellInfo.lateCycRespCells;
    subplot 241
    x = cellInfo.detectComboWeight(ind);
    y = cellInfo.detectWeight{visualTrials}(ind);
    plot(x,y,'.','MarkerSize',10)
    [c,p] = corrcoef(x,y);
    hold on
    figXAxis([],'Combo Weight',weightLim)
    figYAxis([],'Vis. Only Weight',weightLim)
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    title(sprintf('Visual Detect: R=%s,p=%s',...
        num2str(round(c(1,2),2,'significant')),...
        num2str(round(p(1,2),2,'significant'))))
    subplot 242
    x = cellInfo.detectComboWeight(ind);
    y = cellInfo.detectWeight{auditoryTrials}(ind);
    plot(x,y,'.','MarkerSize',10)
    [c,p] = corrcoef(x,y);
    hold on
    figXAxis([],'Combo Weight',weightLim)
    figYAxis([],'Aud. Only Weight',weightLim)
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    title(sprintf('Auditory Detect: R=%s,p=%s',...
        num2str(round(c(1,2),2,'significant')),...
        num2str(round(p(1,2),2,'significant'))))
    subplot 243
    x = cellInfo.detectWeight{visualTrials}(ind);
    y = cellInfo.detectWeight{auditoryTrials}(ind);
    plot(x,y,'.','MarkerSize',10)
    [c,p] = corrcoef(x,y);
    hold on
    figXAxis([],'Vis. Only Weight',weightLim)
    figYAxis([],'Aud. Only Weight',weightLim)
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    title(sprintf('Single Detect: R=%s,p=%s',...
        num2str(round(c(1,2),2,'significant')),...
        num2str(round(p(1,2),2,'significant'))))
    subplot 244
    y = abs(cat(1,cellInfo.detectWeight{visualTrials}(ind),...
        cellInfo.detectWeight{auditoryTrials}(ind),...
        cellInfo.detectComboWeight(ind)));
    yerr = ste(y,2);
    errorbar(1:3,mean(y,2),yerr,'.')
    [p,~,stats] = anova1(y',[],'off');
    c = multcompare(stats,[],'off');
    disp(c(:,[1,2,6]))
    figXAxis([],'Train',[0 4],1:3,{'Visual','Auditory','Vis+Aud'})
    figYAxis([],'Weight Magnitude',[0 0.5])
    figAxForm
    title(sprintf('Detect Model, p=%s',num2str(round(p,2,'significant'))))
    subplot 245
    x = cellInfo.targetComboWeight(ind);
    y = cellInfo.targetWeight{visualTrials}(ind);
    plot(x,y,'.','MarkerSize',10)
    [c,p] = corrcoef(x,y);
    hold on
    figXAxis([],'Combo Weight',weightLim)
    figYAxis([],'Vis. Only Weight',weightLim)
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    title(sprintf('Visual Target: R=%s,p=%s',...
        num2str(round(c(1,2),2,'significant')),...
        num2str(round(p(1,2),2,'significant'))))
    subplot 246
    x = cellInfo.targetComboWeight(ind);
    y = cellInfo.targetWeight{auditoryTrials}(ind);
    plot(x,y,'.','MarkerSize',10)
    [c,p] = corrcoef(x,y);
    hold on
    figXAxis([],'Combo Weight',weightLim)
    figYAxis([],'Aud. Only Weight',weightLim)
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    title(sprintf('Auditory Target: R=%s,p=%s',...
        num2str(round(c(1,2),2,'significant')),...
        num2str(round(p(1,2),2,'significant'))))
    subplot 247
    x = cellInfo.targetWeight{visualTrials}(ind);
    y = cellInfo.targetWeight{auditoryTrials}(ind);
    plot(x,y,'.','MarkerSize',10)
    [c,p] = corrcoef(x,y);
    hold on
    figXAxis([],'Vis. Only Weight',weightLim)
    figYAxis([],'Aud. Only Weight',weightLim)
    figAxForm
    hline(0,'k:')
    vline(0,'k:')
    title(sprintf('Single Target: R=%s,p=%s',...
        num2str(round(c(1,2),2,'significant')),...
        num2str(round(p(1,2),2,'significant'))))
    subplot 248
    y = abs(cat(1,cellInfo.targetWeight{visualTrials}(ind),...
        cellInfo.targetWeight{auditoryTrials}(ind),...
        cellInfo.targetComboWeight(ind)));
    yerr = ste(y,2);
    errorbar(1:3,mean(y,2),yerr,'.')
    [p,~,stats] = anova1(y',[],'off');
    c = multcompare(stats,[],'off');
    disp(c(:,[1,2,6]))
    figXAxis([],'Train',[0 4],1:3,{'Visual','Auditory','Vis+Aud'})
    figYAxis([],'Weight Magnitude',[0 0.5])
    figAxForm
    title(sprintf('Target Model, p=%s',num2str(round(p,2,'significant'))))
    print([fnout 'decodeModelComboTrainWeights'],'-dpdf','-fillpage')

%% stats for modeling
avName = {'Vis','Aud'};
for iav = 1:2
    nav = avName{iav};
    if iav == 2
        stimID = [1,3];
    else
        stimID = 1:3;
    end
    if iav == 1
        [p,~,stats] = anova1(dcModel(1).av(iav).pctCorrect(stimID,:)',[],'off');
        tbl = multcompare(stats,[],'off');
        tbl = tbl(:,[1:2,end]);
    else
        y = dcModel(1).av(iav).pctCorrect(stimID,:)';
        [~,p] = ttest(y(:,1),y(:,2));
        tbl = [];
    end
    imgStats.av(iav).targetPerformanceXStimAnova = p;
    imgStats.av(iav).targetPerformanceXStimPostHocTests = tbl;
    fprintf('%s Target Performance x Stim ANOVA, p=%s\n',nav,num2str(round(...
        imgStats.av(iav).targetPerformanceXStimAnova,2,'significant')))
    disp(imgStats.av(iav).targetPerformanceXStimPostHocTests)
    
    if iav == 1
        [p,~,stats] = anova1(dcModel(2).av(iav).pctCorrect(stimID,:)',[],'off');
        tbl = multcompare(stats,[],'off');
        tbl = tbl(:,[1:2,end]);
    else
        y = dcModel(2).av(iav).pctCorrect(stimID,:)';
        [~,p] = ttest(y(:,1),y(:,2));
        tbl = [];
    end
    imgStats.av(iav).detectPerformanceXStimAnova = p;
    imgStats.av(iav).detectPerformanceXStimPostHocTests = tbl;
    fprintf('%s Detect Performance x Stim ANOVA, p=%s\n',nav,num2str(round(...
        imgStats.av(iav).detectPerformanceXStimAnova,2,'significant')))
    disp(imgStats.av(iav).detectPerformanceXStimPostHocTests)
end

y = [];
z = [];
for iexp = 1:nexp
    if ~isempty(decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget)
        y = cat(1,y,[decodeAnalysis(iexp).av(visualTrials).validMatchedPctCorrectTarget_holdout,...
            decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget]);
    z = cat(1,z,decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectTarget_testAudModel);
    end
end
[~,p] = ttest(y(:,1),y(:,2));
imgStats.visInvVsVisValTargetTest = p;
[~,p] = ttest(z,0.5);
imgStats.visInvTestAudTarget = p;
fprintf('Paired t-Test Vis. Val vs. Vis. Inv. in Vis Target Model,p=%s\n',...
    num2str(round(imgStats.visInvVsVisValTargetTest,2,'significant')))
fprintf('t-Test Vis. Inv. vs. 0.5 in Aud Target Model,p=%s\n',...
    num2str(round(imgStats.visInvTestAudTarget,2,'significant')))

y = [];
z = [];
a = [];
for iexp = 1:nexp
    if ~isempty(decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect)
        y = cat(1,y,[decodeAnalysis(iexp).av(visualTrials).validMatchedPctCorrectDetect_holdout,...
            decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect]);
        z = cat(1,z,decodeAnalysis(iexp).av(visualTrials).invalidPctCorrectDetect_testAudModel);
        a = cat(1,a,decodeAnalysis(iexp).av(visualTrials).validMatchedPctCorrectDetect_testAudModel);
    end
end
[~,p] = ttest(y(:,1),y(:,2));
imgStats.visInvVsVisValDetectTest = p;
[~,p] = ttest(z,y(:,1));
imgStats.visInvTestAudVsVisInvTestVisDetect = p;
[~,p] = ttest(y(:,1),z);
imgStats.visInvTestAudVsVisValTestVisDetect = p;
fprintf('Paired t-Test Vis. Val vs. Vis. Inv. in Vis Detect Model,p=%s\n',...
    num2str(round(imgStats.visInvVsVisValDetectTest,2,'significant')))
fprintf('Paired t-Test Vis. Val. in Vis Detect vs. Vis. Inv. in Aud Detect Model,p=%s\n',...
    num2str(round(imgStats.visInvTestAudVsVisValTestVisDetect,2,'significant')))

save([fnout 'imgStats'],'imgStats')
end
%% Naive decoding figure

if ~strcmp(ds,'FSAV_attentionV1')
    
    dcModel = struct;
    dcModel(1).cellTypeName = {'Tar.';'Dist.';'N.P.'};
    dcModel(1).name = 'Target';
%     cellTypeIDAllCells = sum(cat(2,(cellInfo.lateCycRespCells & ~cellInfo.targetRespCells).*1,...
%         (cellInfo.lateCycRespCells & cellInfo.targetRespCells).*2,...
%         (~cellInfo.lateCycRespCells & cellInfo.targetRespCells).*3),2);
%     dcModel(1).cellTypeID = cellTypeIDAllCells(cellInfo.dcModelCells);
    dcModel.av(1).name = 'Visual';
    dcModel.av(2).name = 'Auditory';
    for iav = 1:2
        dcModel.av(iav).dv = nan(1,nexp);
        dcModel.av(iav).correlation = [];
        dcModel.av(iav).weight = [];
        dcModel.av(iav).pctCorrect = nan(4,nexp);
        dcModel.av(iav).testPerformance = nan(4,nexp);
        dcModel.av(iav).pctCorrect_otherAV = nan(1,nexp);
    end
    for iav = 1:2
        for iexp = 1:nexp
            d = decodeAnalysis(iexp).av(iav);
            dcModel(1).av(iav).dv(iexp) = d.dvTarget;

            dcModel(1).av(iav).correlation = cat(2, ...
                dcModel(1).av(iav).correlation,d.correlationTarget);
            dcModel(1).av(iav).weight = cat(1, ...
                dcModel(1).av(iav).weight,d.weightTarget);
            dcModel(1).av(iav).pctCorrect(:,iexp) = cat(2,...
                d.pctCorrectXStimTarget_holdout,d.pctCorrectAllTarget_holdout);
            dcModel(1).av(iav).testPerformance(:,iexp) = cat(2,...
                d.pctCorrectXStimTarget_train,d.pctCorrectAllTarget_train) - ...
                cat(2,...
                d.pctCorrectXStimTarget_holdout,d.pctCorrectAllTarget_holdout);
            dcModel(1).av(iav).pctCorrect_otherAV(iexp) = d.pctCorrectTarget_otherAV;
        end
    end

    stimLabel = {'0';'HT';'ET';'All'};
    mdlXStimAlpha = 0.05./3;
    
        iexp=1;
    for iav = 1:2
        dc = decodeDataExpt(iexp).av(iav);
        ind = decodeAnalysis(iexp).cellInd;
    %     trialInd = dcTrialsUsed{iexp,iav};
    %     trOut = dc.outcome(trialInd);
        trOut = dc.outcome;

    %     r = dc.tc(:,ind,trialInd);
        r = dc.tc(:,ind,:);
        [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut);
        wTarget = decodeAnalysis(iexp).av(iav).weightTarget;
        [~,targetWeightSort] = sort(wTarget);

        rTarget = {mean(r(:,:,targetTrInd==0),3);mean(r(:,:,targetTrInd==1),3)};
        rTargetErr = {ste(r(:,:,targetTrInd==0),3);ste(r(:,:,targetTrInd==1),3)};

        tfig_d = figure;
        suptitle(sprintf('%s, %s, Dist. Resp, Target Weight Sorted',...
            decodeDataExpt(iexp).exptName,avName{iav}))
        tfig_t = figure;
        suptitle(sprintf('%s, %s, Tar. Resp, Target Weight Sorted',...
            decodeDataExpt(iexp).exptName,avName{iav}))
        for icell = 1:sum(ind)
            figure(tfig_d);
            subplot(5,3,icell)
            cInd = targetWeightSort(icell);
    %         plot(tt_cycTC,rTarget{1}(tcStartFrame:nFr_cyc,cInd),'k-')
            shadedErrorBar_chooseColor(tt_cycTC,rTarget{1}(tcStartFrame:nFr_cyc,cInd),...
                rTargetErr{1}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
            figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
            figYAxis([],'dF/F',cellRespTCLim)  
            hline(0,'r:')
            vline(respWinTT,'k--')
            figAxForm
            title(sprintf('D Resp., Cell #%s, W=%s',num2str(cInd),...
                num2str(round(wTarget(cInd),2,'significant'))))
            figure(tfig_t);
            subplot(5,3,icell)
    %         plot(tt_cycTC,rTarget{2}(tcStartFrame:nFr_cyc,cInd),'k-')
            shadedErrorBar_chooseColor(tt_cycTC,rTarget{2}(tcStartFrame:nFr_cyc,cInd),...
                rTargetErr{2}(tcStartFrame:nFr_cyc,cInd),[0 0 0]);
            figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],[0,200],[0,200])
            figYAxis([],'dF/F',cellRespTCLim)  
            hline(0,'r:')
            vline(respWinTT,'k--')
            figAxForm
            title(sprintf('T Resp., Cell #%s, W=%s',num2str(cInd),...
                num2str(round(wTarget(cInd),2,'significant'))))
        end
        figure(tfig_d);
        print(...
            [fnout 'targetModel_' avName{iav}(1:3) '_distResp_' decodeDataExpt(iexp).exptName],...
            '-dpdf','-fillpage')
        figure(tfig_t);
        print(...
            [fnout 'targetModel_' avName{iav}(1:3) '_tarResp_' decodeDataExpt(iexp).exptName],...
            '-dpdf','-fillpage')
    end

    
    figure
    for iav = 1:2
        y = dcModel.av(iav).pctCorrect;
        yerr = ste(dcModel.av(iav).pctCorrect,2);
        stimInd = sum(~isnan(y),2) > minTrN;
        y = y(stimInd,:);
%         yerr = yerr(stimInd,:);
        x = 1:sum(stimInd);
        [mdlTest,mdlP] = ttest(y,0.5,'dim',2,'tail','right','alpha',mdlXStimAlpha);
        subplot(1,2,iav)
        for iexp = 1:nexp
            ind = ~isnan(y(:,iexp));
            hold on
            plot(x(ind),y(ind,iexp),'k.-');
        end
        for itest = 1:length(x)
            if mdlTest(itest) == 1
                text(x(itest),1,sprintf('*%s',...
                    num2str(round(mdlP(itest),2,'significant'))))
            else
                text(x(itest),1,sprintf('%s',...
                    num2str(round(mdlP(itest),2,'significant'))))
            end
        end
        errorbar(x,nanmean(y,2),yerr,'.-','MarkerSize',15)
        figXAxis([],'Stim',[0 length(x)+1],1:length(x),stimLabel(stimInd))
        hline(0.5,'k:')
        figYAxis([],'Train - HO Frac. Correct',[0 1])
        figAxForm
        title(sprintf('%s Model, Train %s, Test %s',dcModel.name,...
            dcModel.av(iav).name,dcModel.av(iav).name))
    end
    
    print([fnout 'decodeModelPctCorrectXStim'],'-dpdf','-fillpage')
    
    figure
    ind = cellInfo.lateCycRespCells & ~isnan(cellInfo.targetWeight{visualTrials})';
    subplot 221
    x = antiAnalysis.lateCycSI(ind)';
    y = cellInfo.targetWeight{visualTrials}(ind)';
    hold on
    plot(x,y,'.')
    [corrmat,pmat] = corrcoef(x,y);
    c = corrmat(1,2);
    p = pmat(1,2);
    mdl = fitlm(x,y);
    rsq = mdl.Rsquared.Ordinary;
    yfit = predict(mdl,x);
    plot(x,yfit,'-')
    figXAxis([],'Selectivity',siLim)
    figYAxis([],'Visual Weight',weightLim)
    figAxForm
    title(sprintf('Visual Target Model, Rsq = %s, corr = %s, p = %s',num2str(round(rsq,2,'significant')),...
        num2str(round(c,2,'significant')),num2str(round(p,2,'significant'))))
    subplot 222
    y = cellInfo.targetWeight{auditoryTrials}(ind)';
    hold on
    plot(x,y,'.')
    [corrmat,pmat] = corrcoef(x,y);
    c = corrmat(1,2);
    p = pmat(1,2);
    mdl = fitlm(x,y);
    rsq = mdl.Rsquared.Ordinary;
    yfit = predict(mdl,x);
    plot(x,yfit,'-')
    figXAxis([],'Selectivity',siLim)
    figYAxis([],'Auditory Weight',weightLim)
    figAxForm
    title(sprintf('Auditory Target Model, Rsq = %s, corr = %s, p = %s',num2str(round(rsq,2,'significant')),...
        num2str(round(c,2,'significant')),num2str(round(p,2,'significant'))))
    subplot 223
    ind = cellInfo.lateCycRespCells & ~isnan(cellInfo.targetWeight{visualTrials})';
    x = cellInfo.targetWeight{visualTrials}(ind)';
    y = cellInfo.targetWeight{auditoryTrials}(ind)';
    hold on
    plot(x,y,'.')
    [corrmat,pmat] = corrcoef(x,y);
    c = corrmat(1,2);
    p = pmat(1,2);
    mdl = fitlm(x,y);
    rsq = mdl.Rsquared.Ordinary;
    yfit = predict(mdl,x);
    plot(x,yfit,'-')
    figXAxis([],'Visual Weight',weightLim)
    figYAxis([],'Auditory Weight',weightLim)
    figAxForm
    title(sprintf('Target Model, Rsq = %s, corr = %s, p = %s',num2str(round(rsq,2,'significant')),...
        num2str(round(c,2,'significant')),num2str(round(p,2,'significant'))))
    
    print([fnout 'decodeModelWeights'],'-dpdf','-fillpage')
end