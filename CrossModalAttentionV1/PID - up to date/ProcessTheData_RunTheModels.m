% cd('...\For Jeff') % fill this in
load('i613_i614_i625_i668_trOutcomeStruct_cells_attentionV1')
load('attentionV1_200424_imgAnalysisData.mat')
doLoadPreviousAnalysis = false;
doDecoding = true;
%%
imgParams_FSAV
bxParams_FSAV_attnV1ms

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

trOutType = {'h';'m';'fa';'cr';'yes';'no'};
% trOutTypeName = {'H-All';'H-HT';'H-ET';'M-All';'M-HT';'M-ET';'FA';'CR'};
% trOutTypeName = {'H';'M';'FA';'CR'};
%% pool experiment data
rng(0)
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
        if imouse == 1 & iexp == 1
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

        maxCycles = min([max(d.av(visualTrials).align(alignStart).nCycles),...
            max(d.av(auditoryTrials).align(alignStart).nCycles)]);

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
                elseif iav == 1 | strcmp(ds,'FSAV_V1_audControl')
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
            trOutTC = cell(1,length(trOutType));

            if strcmp(ds,'FSAV_attentionV1') | strcmp(ds,'FSAV_attentionV1_noAttn')
                hitInd = strcmp(trOut,'h');
                missInd = strcmp(trOut,'m');
                matchInd = getMatchedOutcomeTrialIndex(...
                    trOut,trStimID,minTrN);
                trOutTC = {trTC(:,:,matchInd & hitInd),trTC(:,:,matchInd & missInd),...
                    trTC(:,:,strcmp(trOut,'fa')),trTC(:,:,strcmp(trOut,'cr')),...
                    cat(3,trTC(:,:,matchInd & hitInd),trTC(:,:,strcmp(trOut,'fa'))),...
                    cat(3,trTC(:,:,matchInd & missInd),trTC(:,:,strcmp(trOut,'cr')))};
            end
            targetDataExpt(exptN).av(iav).trOutTC = trOutTC;

            de = d.av(iav).align(alignStart);
            if strcmp(ds,'FSAV_attentionV1')|strcmp(ds,'FSAV_V1_audControl')
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
            if iav == 1 | strcmp(ds,'FSAV_V1_audControl')
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
        antiDataExpt(exptN).mouse = mouse(imouse).mouse_name;
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
        eaCycSI = cellfun(@(x,y) ...
            getSelectivityIndex(squeeze(mean(x(respwin,:,:),1)-mean(x(basewin_0,:,:),1))',...
            squeeze(mean(y(respwin,:,:),1)-mean(y(basewin_0,:,:),1))'),...
            cycTC(visualTrials,:),...
            cycTC(auditoryTrials,:),'unif',0);
        antiDataExpt(exptN).eaCycSI = eaCycSI;

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
    respCellsExpt(iexp).mouse = antiDataExpt(iexp).mouse;
    firstTCAV = cat(3,antiDataExpt(iexp).cycTC{visualTrials,1},...
        antiDataExpt(iexp).cycTC{auditoryTrials,1});
    fprintf('%s: %s first stim\n',antiDataExpt(iexp).exptName,num2str(size(firstTCAV,3)))

    lateTC_vis = [];
    lateTC_aud = [];
    lateCycles = lateCyclesMin:size(antiDataExpt(iexp).cycTC,2);
    for icyc = 1:length(lateCycles)
        lateTC_vis = cat(3,lateTC_vis,...
            antiDataExpt(iexp).cycTC{visualTrials,lateCycles(icyc)});
        lateTC_aud = cat(3,lateTC_aud,...
            antiDataExpt(iexp).cycTC{auditoryTrials,lateCycles(icyc)});
    end
    lateTCAV = cat(3,lateTC_vis,lateTC_aud);
    fprintf('%s: %s late stim\n',antiDataExpt(iexp).exptName,num2str(size(lateTCAV,3)))

    longTCAV = cat(3,antiDataExpt(iexp).longTC{visualTrials,1},...
        antiDataExpt(iexp).longTC{auditoryTrials,1});
    fprintf('%s: %s long trials\n',antiDataExpt(iexp).exptName,num2str(size(longTCAV,3)))

    % noise correlation

    rv = squeeze(mean(lateTC_vis(respwin,:,:),1)-mean(lateTC_vis(basewin_0,:,:),1))';
    ra = squeeze(mean(lateTC_aud(respwin,:,:),1)-mean(lateTC_aud(basewin_0,:,:),1))';
    rav = cat(1,rv,ra);
%         rv_rect = rv;
%         rv_rect(rv<0) = nan;
%         ra_rect = ra;
%         ra_rect(ra<0) = nan;
%         rav_rect = cat(1,rv_rect,ra_rect);
%         rv_rect = rv-min(rv)+1;
%         ra_rect = ra-min(ra)+1;
%         rav_rect = cat(1,rv,ra)-min(cat(1,rv,ra))+1;

    rsc_vis = corrcoef(rv);
    rsc_aud = corrcoef(ra);
    rsc_all = corrcoef(rav);

    nc = size(rv,2);
    gm_vis = nan(nc);
    gm_aud = nan(nc);
    gm_all = nan(nc);
    for icell1 = 1:nc
        rv1 = rv(:,icell1);
        ra1 = ra(:,icell1);
        rav1 = rav(:,icell1);
        for icell2 = 1:nc
            rv2 = rv(:,icell2);
            ra2 = ra(:,icell2);
            rav2 = rav(:,icell2);
            d = [rv1;rv2];
            gm_vis(icell1,icell2) = geomean(d-min(d)+1)+min(d)-1;
            d = [ra1;ra2];
            gm_aud(icell1,icell2) = geomean(d-min(d)+1)+min(d)-1;
            d = [rav1;rav2];
            gm_all(icell1,icell2) = geomean(d-min(d)+1)+min(d)-1;
        end
    end

    respCellsExpt(iexp).gm = {gm_vis,gm_aud,gm_all};
    respCellsExpt(iexp).rsc = {rsc_vis,rsc_aud,rsc_all};
    % randomly select trials for responsive cells test so that each
    % experiment test the same number of trials
    rng(0)
    if size(lateTCAV,3) > minTrN_lateresp
        ind = randsample(size(firstTCAV,3),minTrN_firstresp);
        firstRespCells = ttest(...
            squeeze(mean(firstTCAV(respwin,:,ind),1)),...
            squeeze(mean(firstTCAV(basewin_0,:,ind),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);
        ind = randsample(size(lateTCAV,3),minTrN_lateresp);
        lateCycRespCells = ttest(...
            squeeze(mean(lateTCAV(respwin,:,ind),1)),...
            squeeze(mean(lateTCAV(basewin_0,:,ind),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);
    else
        firstRespCells = ttest(...
            squeeze(mean(firstTCAV(respwin,:,:),1)),...
            squeeze(mean(firstTCAV(basewin_0,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);
        lateCycRespCells = ttest(...
            squeeze(mean(lateTCAV(respwin,:,:),1)),...
            squeeze(mean(lateTCAV(basewin_0,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);
    end        
    if size(longTCAV,3) > minTrN_latewin
        ind = randsample(size(longTCAV,3),minTrN_latewin);
        lateRespCells = ttest(...
            squeeze(mean(longTCAV(lateWinFr,:,ind),1)),...
            squeeze(mean(longTCAV(basewin,:,ind),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);
        lateSuppCells = ttest(...
            squeeze(mean(longTCAV(lateWinFr,:,ind),1)),...
            squeeze(mean(longTCAV(basewin,:,ind),1)),...
            'dim',2,'tail','left','alpha',cellGroupsAlpha);
    else
        lateRespCells = ttest(...
            squeeze(mean(longTCAV(lateWinFr,:,:),1)),...
            squeeze(mean(longTCAV(basewin,:,:),1)),...
            'dim',2,'tail','right','alpha',cellGroupsAlpha);
        lateSuppCells = ttest(...
            squeeze(mean(longTCAV(lateWinFr,:,:),1)),...
            squeeze(mean(longTCAV(basewin,:,:),1)),...
            'dim',2,'tail','left','alpha',cellGroupsAlpha);
    end

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
%% do the decoding
if doDecoding
    decodeAnalysis = struct;
    decodeAnalysis.av(visualTrials).name = 'Visual';
    decodeAnalysis.av(auditoryTrials).name = 'Auditory';
    for iexp = 1:nexp
            rng(0) % cells randomized and trials randomized
%             cellInd = respCellsExpt(iexp).decodeAnalysisCells & ...
%                 oriTuningExpt(iexp).tuningReliability' <= tuningReliabilityThresh_decode;
            cellInd = respCellsExpt(iexp).lateCycRespCells|...
                respCellsExpt(iexp).targetRespCells;
            decodeAnalysis(iexp).nCells = sum(cellInd);            
            
%             decodeAnalysis(iexp).nCellsSelected = sum(cellInd);
            decodeAnalysis(iexp).cellInd = cellInd;
            respOther = cell(1,2);
            trOutOther = cell(1,2);
            trSampleIndOther = cell(1,2);
            detectGLMOther = cell(1,2);
            targetGLMOther = cell(1,2);
            if isfield(decodeDataExpt(iexp).av(visualTrials),'catchOutcome')
                data = cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).resp(cellInd,:),...
                    decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,:),...
                    decodeDataExpt(iexp).av(visualTrials).catchResp(cellInd,:))';               
                
                [coeffAllCells,scoresAllCells,latentAllCells] = pca(data);                
                respAllCells_av_withInv = zscore(scoresAllCells);
                nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
                naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
                ninv = size(decodeDataExpt(iexp).av(visualTrials).catchResp,2);
                respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
                respAllCells_aud = respAllCells_av_withInv((nvis+1):(nvis+naud),:);
                catchResp = respAllCells_av_withInv((end-ninv+1):end,:);
                
                data4Shuff = cell(1,2);
                data4Shuff{visualTrials} = data(1:nvis,:);
                data4Shuff{auditoryTrials} = data((nvis+1):(nvis+naud),:);
                                
                stim_eaTr = cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).stim,...
                    decodeDataExpt(iexp).av(auditoryTrials).stim,...
                    decodeDataExpt(iexp).av(visualTrials).catchStim(...
                    1:length(decodeDataExpt(iexp).av(visualTrials).catchOutcome)));
                stims_all = unique(stim_eaTr);
                mod_eaTr = cat(2,...
                    ones(1,length(decodeDataExpt(iexp).av(visualTrials).stim)),...
                    2.*ones(1,length(decodeDataExpt(iexp).av(auditoryTrials).stim)),...
                    3.*ones(1,length(decodeDataExpt(iexp).av(visualTrials).catchOutcome)));
                data_resid = nan(size(data));
                for istim = 1:length(stims_all)
                    for imod = 1:3
                        ind = stim_eaTr == stims_all(istim) & mod_eaTr == imod;
                        data_resid(ind,:) = data(ind,:) - mean(data(ind,:));
                    end
                end
                
                [coeffAllCells_resid,scoresAllCells_resid] = pca(data_resid);                
                respAllCells_av_withInv_resid = zscore(scoresAllCells_resid);
                respAllCells_vis_resid = respAllCells_av_withInv_resid(1:nvis,:);
                respAllCells_aud_resid = respAllCells_av_withInv_resid((nvis+1):(nvis+naud),:);
                catchResp_resid = respAllCells_av_withInv_resid((end-ninv+1):end,:);    
%                 [coeffAllCells_vis,scoresAllCells_vis,latentAllCells_vis] = pca(respAllCells_vis(:,cellInd));
%                 [coeffAllCells_aud,scoresAllCells_aud,latentAllCells_aud] = pca(respAllCells_aud(:,cellInd));
%                 [coeffAllCells_catch,scoresAllCells_catch,latentAllCells_catch] = pca(catchResp(:,cellInd));

            else
                data = cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).resp(cellInd,:),...
                    decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,:))';
                
                [coeffAllCells,scoresAllCells,latentAllCells] = pca(data);
                respAllCells_av_withInv = zscore(scoresAllCells);
                nvis = size(decodeDataExpt(iexp).av(visualTrials).resp,2);
                naud = size(decodeDataExpt(iexp).av(auditoryTrials).resp,2);
                respAllCells_vis = respAllCells_av_withInv(1:nvis,:);
                respAllCells_aud = respAllCells_av_withInv((nvis+1):end,:);

                
                data4Shuff = cell(1,2);
                data4Shuff{visualTrials} = data(1:nvis,:);
                data4Shuff{auditoryTrials} = data((nvis+1):(nvis+naud),:);
%                 [coeffAllCells_vis,scoresAllCells_vis,latentAllCells_vis] = pca(respAllCells_vis(:,cellInd));
%                 [coeffAllCells_aud,scoresAllCells_aud,latentAllCells_aud] = pca(respAllCells_aud(:,cellInd));

                stim_eaTr = cat(2,...
                    decodeDataExpt(iexp).av(visualTrials).stim,...
                    decodeDataExpt(iexp).av(auditoryTrials).stim);
                stims_all = unique(stim_eaTr);
                mod_eaTr = cat(2,...
                    ones(1,length(decodeDataExpt(iexp).av(visualTrials).stim)),...
                    2.*ones(1,length(decodeDataExpt(iexp).av(auditoryTrials).stim)));
                data_resid = nan(size(data));
                for istim = 1:length(stims_all)
                    for imod = 1:2
                        ind = stim_eaTr == stims_all(istim) & mod_eaTr == imod;
                        data_resid(ind,:) = data(ind,:) - mean(data(ind,:),1);
                    end
                end
                
                [coeffAllCells_resid,scoresAllCells_resid] = pca(data_resid);                
                respAllCells_av_withInv_resid = zscore(scoresAllCells_resid);
                respAllCells_vis_resid = respAllCells_av_withInv_resid(1:nvis,:);
                respAllCells_aud_resid = respAllCells_av_withInv_resid((nvis+1):(nvis+naud),:); 
            end
            rng(0)
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;

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
                        respStimSort{istim} = respAllCells(indSample,:);
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
                decodeAnalysis(iexp).av(iav).respAllCells = respAllCells(matchTrialsInd,:);
                decodeAnalysis(iexp).av(iav).trOut = trOut(matchTrialsInd);
                
                respStimSort_av{iav} = respStimSort;
                trOutStimSort_av{iav} = trOutStimSort;
                stimSortInd_av{iav} = stimSortInd;
                matchTrialsInd_av{iav} = matchTrialsInd;
            end
            
            data4Shuff_matched = cellfun(@(x,y) x(y,:),data4Shuff,matchTrialsInd_av,'unif',0);
            trOut_matched = cellfun(@(x,y) x(y),...
                {decodeDataExpt(iexp).av(visualTrials).outcome,...
                decodeDataExpt(iexp).av(auditoryTrials).outcome},matchTrialsInd_av,'unif',0);
            [detectTrInd,targetTrInd] = cellfun(@(x) getStimAndBehaviorYs(x),...
                trOut_matched,'unif',0);
            dataShuff_detect = cellfun(@(x,y) shuffWithinCats(x,y),data4Shuff_matched,...
                detectTrInd,'unif',0);
            dataShuff_target = cellfun(@(x,y) shuffWithinCats(x,y),data4Shuff_matched,...
                targetTrInd,'unif',0);
            
            n = cellfun(@(x) size(x,1),data4Shuff_matched);
            d = cell2mat(cellfun(@(x) cell2mat(x'),dataShuff_detect,'unif',0)');
            [~,scoresAllCells] = pca(d);
            pcaShuff_zscore = zscore(scoresAllCells);
            pcaShuff_detect_matched = {pcaShuff_zscore(1:n(1),:),pcaShuff_zscore((n(1)+1):end,:)};
            d = cell2mat(cellfun(@(x) cell2mat(x'),dataShuff_target,'unif',0)');
            [~,scoresAllCells] = pca(d);
            pcaShuff_zscore = zscore(scoresAllCells);
            pcaShuff_target_matched = {pcaShuff_zscore(1:n(1),:),pcaShuff_zscore((n(1)+1):end,:)};
            
            nTrials = min(cellfun(@length,matchTrialsInd_av));
            if round(nTrials.*fracPCsUsed) > maxPCs
                if round(nTrials.*fracPCsUsed) > sum(cellInd)
                    nPCs = sum(cellInd);
                else
                    nPCs = maxPCs;
                end
            elseif sum(cellInd) < round(nTrials.*fracPCsUsed)
                nPCs = sum(cellInd);
            else
                nPCs = round(nTrials.*fracPCsUsed);
            end
            decodeAnalysis(iexp).nPCs = nPCs;
            decodeAnalysis(iexp).minTrials = nTrials;
            
            for iav = 1:2
                trOut = decodeDataExpt(iexp).av(iav).outcome;
                respStimSort = respStimSort_av{iav};
                trOutStimSort = trOutStimSort_av{iav};
                stimSortInd = stimSortInd_av{iav};
                matchTrialsInd = matchTrialsInd_av{iav};

                if iav == 1
                    respAllCells = respAllCells_vis;
%                     trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                elseif iav == 2
                    respAllCells = respAllCells_aud;
%                     trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                end
                
                emptyInd = cellfun(@isempty,respStimSort);
                respStimSort(~emptyInd) = cellfun(@(x) x(:,1:nPCs),...
                    respStimSort(~emptyInd),'unif',0);
                
                resp = respAllCells(matchTrialsInd,1:nPCs);
%                 resp = respAllCells(matchTrialsInd,cellInd);
                [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));

                detectCorr = corr(detectTrInd,resp);
                targetCorr = corr(targetTrInd,resp);

                C = eye(size(resp,2));
                p=1;
                [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
                [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');
                                
                detectWeight = detectGLM.beta(2:end);
                targetWeight = targetGLM.beta(2:end);

                dv_detect = mean(detectTrInd);
                dv_target = mean(targetTrInd);
                
                fprintf('Expt %s, hold-out analysis\n',num2str(iexp))
                pctCorrectDetect_train = getPctCorr_trainData(detectGLM,resp,detectTrInd,dv_detect);
                pctCorrectDetect_ho = getPctCorr_hoData(resp,detectTrInd,dv_detect);

                pctCorrectTarget_train = getPctCorr_trainData(targetGLM,resp,targetTrInd,dv_target);
                pctCorrectTarget_ho = getPctCorr_hoData(resp,targetTrInd,dv_target);
 
                pctCorrDetect_xStim_train = nan(1,nStimBins);
                pctCorrDetect_xStim_ho = nan(1,nStimBins);
                pctCorrTarget_xStim_train = nan(1,nStimBins);
                pctCorrTarget_xStim_ho = nan(1,nStimBins);
                for istim = 1:nStimBins
                        if isempty(trOutStimSort{istim})
                            continue
                        end
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

                fprintf('Expt %s, starting resp win analysis\n',num2str(iexp))
                nwins = size(decodeDataExpt(iexp).av(iav).movWinResp,3);
                pctCorrectDetect_train_respwin = nan(1,nwins);
                pctCorrectDetect_ho_respwin = nan(1,nwins);
                pctCorrectTarget_train_respwin = nan(1,nwins);
                pctCorrectTarget_ho_respwin = nan(1,nwins);
                for iwin = 1:nwins
%                     rAllCells = zscore(decodeDataExpt(iexp).av(iav).movWinResp(:,:,iwin)');
%                     r = rAllCells(matchTrialsInd,cellInd);
                    rAllCells = zscore(decodeDataExpt(iexp).av(iav).movWinResp(:,:,iwin)');
                    [~,sAllCells] = pca(rAllCells(:,cellInd));
                    r = sAllCells(matchTrialsInd,1:nPCs);
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
                
                % train the model on the PCs from the residuals
                if iav == 1
                    respAllCells_resid = respAllCells_vis_resid;
                elseif iav == 2
                    respAllCells_resid = respAllCells_aud_resid;
                end
                resp_resid = respAllCells_resid(matchTrialsInd,1:nPCs);
                C = eye(size(resp_resid,2));
                [~,~,detectGLM_resid] = glmfit(resp_resid*C,detectTrInd,'binomial');
                [~,~,targetGLM_resid] = glmfit(resp_resid*C,targetTrInd,'binomial');
                 
                pctCorrectDetect_train_resid = getPctCorr_trainData(detectGLM_resid,resp_resid,detectTrInd,dv_detect);
                pctCorrectDetect_ho_resid = getPctCorr_hoData(resp_resid,detectTrInd,dv_detect);

                pctCorrectTarget_train_resid = getPctCorr_trainData(targetGLM_resid,resp_resid,targetTrInd,dv_target);
                pctCorrectTarget_ho_resid = getPctCorr_hoData(resp_resid,targetTrInd,dv_target);
 
                detectWeight_cells_resid = coeffAllCells_resid(:,1:nPCs)*detectGLM_resid.beta(2:(nPCs+1));
                targetWeight_cells_resid = coeffAllCells_resid(:,1:nPCs)*targetGLM_resid.beta(2:(nPCs+1));
                
                % train the model with stim or choices as predictors
                resp_withStim = cat(2,resp,targetTrInd);
                resp_withChoice = cat(2,resp,detectTrInd);
                C = eye(size(resp_withStim,2));
                [~,~,detectGLM_withStim] = glmfit(resp_withStim*C,detectTrInd,'binomial');
                pctCorrectDetect_ho_withStim = getPctCorr_hoData(resp_withStim,detectTrInd,dv_detect);
                C = eye(size(resp_withChoice,2));
                [~,~,targetGLM_withChoice] = glmfit(resp_withChoice*C,targetTrInd,'binomial');
                pctCorrectTarget_ho_withChoice = getPctCorr_hoData(resp_withChoice,targetTrInd,dv_target);
                
                detectWeight_cells_withStim = coeffAllCells(:,1:nPCs)*detectGLM_withStim.beta(2:(nPCs+1));
                targetWeight_cells_withChoice = coeffAllCells(:,1:nPCs)*targetGLM_withChoice.beta(2:(nPCs+1));
                modeledStimWeight = detectGLM_withStim.beta(nPCs+2);
                modeledChoiceWeight = targetGLM_withChoice.beta(nPCs+2);
                
%                 pctCorrectDetect_ho_withStim_pcsOnly = getPctCorr_hoData_choosePCs(...
%                     resp_withStim,detectTrInd,dv_detect,1:nPCs);
%                 pctCorrectTarget_ho_withChoice_pcsOnly = getPctCorr_hoData_choosePCs(...
%                     resp_withChoice,targetTrInd,dv_target,1:nPCs);
                
                
                % train the model with fixed stim or choice predictors
                weight_choice = glmfit(detectTrInd,targetTrInd,'binomial');
                weight_stim = glmfit(targetTrInd,detectTrInd,'binomial');
                
                dv_choiceEaTrial = glmval(weight_choice,detectTrInd,'logit');
                dv_stimEaTrial = glmval(weight_stim,targetTrInd,'logit');
                
                resp_withFixedChoice = cat(2,resp,dv_choiceEaTrial);
                resp_withFixedStim = cat(2,resp,dv_stimEaTrial);
                
                C = eye(size(resp_withFixedChoice,2));
                [~,~,targetGLM_withFixedChoice] = glmfit(resp_withFixedChoice*C,targetTrInd,'binomial');
                pctCorrectTarget_ho_withFixedChoice = getPctCorr_hoData(resp_withFixedChoice,targetTrInd,dv_target);
                C = eye(size(resp_withFixedStim,2));
                [~,~,detectGLM_withFixedStim] = glmfit(resp_withFixedStim*C,detectTrInd,'binomial');
                pctCorrectDetect_ho_withFixedStim = getPctCorr_hoData(resp_withFixedStim,detectTrInd,dv_detect);
                
                detectGLMOther_withStim{iav} = detectGLM_withFixedStim;
                targetGLMOther_withChoice{iav} = targetGLM_withFixedChoice; 
                
                pctCorrectDetect_ho_withFixedStim_pcsOnly = getPctCorr_hoData_choosePCs(...
                    resp_withFixedStim,detectTrInd,dv_detect,1:nPCs);
                pctCorrectTarget_ho_withFixedChoice_pcsOnly = getPctCorr_hoData_choosePCs(...
                    resp_withFixedChoice,targetTrInd,dv_target,1:nPCs);
                
                
                
%                 C = eye(size(dv_stimEaTrial,2));
%                 [~,~,detectGLM_stimOnly] = glmfit(dv_stimEaTrial*C,detectTrInd,'binomial');
                pctCorrectDetect_stimOnly = getPctCorr_hoData(dv_stimEaTrial,detectTrInd,dv_detect);
%                 C = eye(size(dv_choiceEaTrial,2));
%                 [~,~,targetGLM_choiceOnly] = glmfit(dv_choiceEaTrial*C,detectTrInd,'binomial');
                pctCorrectTarget_choiceOnly = getPctCorr_hoData(dv_choiceEaTrial,targetTrInd,dv_target);
                
%                 % bootstrap trials used to generate stim:choice correlation
%                 % distribution
%                 choiceModel_fixedStim_struct = struct;
%                 [choiceModel_fixedStim_struct.stimChoiceCorr_boot,...
%                     choiceModel_fixedStim_struct.trialsUsed_boot,...
%                     choiceModel_fixedStim_struct.choiceModelWeights_boot,...
%                     choiceModel_fixedStim_struct.pctCorr_choice_fixedStim_boot,...
%                     choiceModel_fixedStim_struct.pctCorr_choice_fixedStim_pcsOnly_boot,...
%                     choiceModel_fixedStim_struct.pctCorr_choice_fixedStimOnly_boot,...
%                     choiceModel_fixedStim_struct.pctCorr_pcs,...
%                     choiceModel_fixedStim_struct.weights_pcs] = ...
%                     getSubSampledModel_boot(...
%                     targetTrInd,detectTrInd,resp_withFixedStim,100,round(0.5*length(detectTrInd)));
%                 choiceModel_fixedStim_struct.choiceModelWeights_boot = ...
%                     cellfun(@(x,y) cat(1,coeffAllCells(:,1:nPCs)*x(2:(end-1)),x(end)),...
%                     choiceModel_fixedStim_struct.choiceModelWeights_boot,'unif',0);
%                 choiceModel_fixedStim_struct.weights_pcs = ...
%                     cellfun(@(x,y) coeffAllCells(:,1:nPCs)*x(2:(end)),...
%                     choiceModel_fixedStim_struct.weights_pcs,'unif',0);
%                 
%                 stimModel_fixedChoice_struct = struct;
%                 [stimModel_fixedChoice_struct.stimChoiceCorr_boot,...
%                     stimModel_fixedChoice_struct.trialsUsed_boot,...
%                     stimModel_fixedChoice_struct.stimModelWeights_boot,...
%                     stimModel_fixedChoice_struct.pctCorr_stim_fixedChoice_boot,...
%                     stimModel_fixedChoice_struct.pctCorr_stim_fixedChoice_pcsOnly_boot,...
%                     stimModel_fixedChoice_struct.pctCorr_stim_fixedChoiceOnly_boot,...
%                     stimModel_fixedChoice_struct.pctCorr_pcs,...
%                     stimModel_fixedChoice_struct.weights_pcs] = ...
%                     getSubSampledModel_boot(...
%                     detectTrInd,targetTrInd,resp_withFixedChoice,100,round(0.5*length(targetTrInd)));
%                 stimModel_fixedChoice_struct.stimModelWeights_boot = ...
%                     cellfun(@(x,y) cat(1,coeffAllCells(:,1:nPCs)*x(2:(end-1)),x(end)),...
%                     stimModel_fixedChoice_struct.stimModelWeights_boot,'unif',0);
%                 stimModel_fixedChoice_struct.weights_pcs = ...
%                     cellfun(@(x,y) coeffAllCells(:,1:nPCs)*x(2:(end)),...
%                     stimModel_fixedChoice_struct.weights_pcs,'unif',0);
                
                % train trial shuffle model
                detectTrInd_shuff = cat(1,zeros(length(detectTrInd)-sum(detectTrInd),1),...
                    ones(sum(detectTrInd),1));
                targetTrInd_shuff = cat(1,zeros(length(targetTrInd)-sum(targetTrInd),1),...
                    ones(sum(targetTrInd),1));
                resp_detect_shuff = pcaShuff_detect_matched{iav}(:,1:nPCs);
                resp_target_shuff = pcaShuff_target_matched{iav}(:,1:nPCs);
                
                detectCorr_shuff = corr(detectTrInd_shuff,resp_detect_shuff);
                targetCorr_shuff = corr(targetTrInd_shuff,resp_target_shuff);

                C = eye(size(resp_detect_shuff,2));
                p=1;
                [~,~,detectGLM_shuff] = glmfit(resp_detect_shuff*C,detectTrInd,'binomial');
                [~,~,targetGLM_shuff] = glmfit(resp_target_shuff*C,targetTrInd,'binomial');
                                
                detectWeight_shuff = detectGLM_shuff.beta(2:end);
                targetWeight_shuff = targetGLM_shuff.beta(2:end);

                dv_detect_shuff = mean(detectTrInd_shuff);
                dv_target_shuff = mean(targetTrInd_shuff);
                
                fprintf('Expt %s, trial shuffle analysis\n',num2str(iexp))
                pctCorrectDetect_train_shuff = getPctCorr_trainData(detectGLM_shuff,...
                    resp_detect_shuff,detectTrInd_shuff,dv_detect_shuff);
                pctCorrectDetect_ho_shuff = getPctCorr_hoData(resp_detect_shuff,...
                    detectTrInd_shuff,dv_detect_shuff);
                pctCorrectTarget_train_shuff = getPctCorr_trainData(targetGLM_shuff,...
                    resp_target_shuff,targetTrInd_shuff,dv_target_shuff);
                pctCorrectTarget_ho_shuff = getPctCorr_hoData(resp_target_shuff,...
                    targetTrInd_shuff,dv_target_shuff);
                                
                respOther_shuff{iav} = {resp_detect_shuff,resp_target_shuff};
                trOutInd_shuff{iav} = {detectTrInd_shuff,targetTrInd_shuff};
                detectGLMOther_shuff{iav} = detectGLM_shuff;
                targetGLMOther_shuff{iav} = targetGLM_shuff;  
                
                
                % train detect model with distractor or target choices only
                if iav == 1
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                elseif iav == 2
                    trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,ampBins);
                end
%                 trOut_matched = trOut(matchTrialsInd);
                % distractors
                distInd = strcmp(trOut,'fa')|strcmp(trOut,'cr');
                trOut_dist = trOut(distInd);
                resp_dist = respAllCells(distInd,1:nPCs);
                [detectTrInd_dist, ~] = getStimAndBehaviorYs(trOut_dist);
                dv_detect_dist = mean(detectTrInd_dist);
                C = eye(size(resp_dist,2));
                p=1;
                [~,~,detectGLM_dist] = glmfit(resp_dist*C,detectTrInd_dist,'binomial');
                pctCorrectDetect_ho_dist = getPctCorr_hoData(resp_dist,detectTrInd_dist,dv_detect_dist);
                respOther_dist{iav} = resp_dist;
                detectGLMOther_dist{iav} = detectGLM_dist;
                trOutOther_dist{iav} = trOut_dist;
                detectWeight_cells_distOnly = coeffAllCells(:,1:nPCs)*detectGLM_dist.beta(2:end);
                
                % targets
                tarInd = strcmp(trOut,'h')|strcmp(trOut,'m');
                trOut_tar = trOut(tarInd);
                resp_tar = respAllCells(tarInd,1:nPCs);
                
                [detectTrInd_tar, ~] = getStimAndBehaviorYs(trOut_tar);
                dv_detect_tar = mean(detectTrInd_tar);
                C = eye(size(resp_tar,2));
                p=1;
                [~,~,detectGLM_tar] = glmfit(resp_tar*C,detectTrInd_tar,'binomial');
                pctCorrectDetect_ho_tar = getPctCorr_hoData(resp_tar,detectTrInd_tar,dv_detect_tar);
                respOther_tar{iav} = resp_tar;
                detectGLMOther_tar{iav} = detectGLM_tar;
                trOutOther_tar{iav} = trOut_tar;
                detectWeight_cells_tarOnly = coeffAllCells(:,1:nPCs)*detectGLM_tar.beta(2:end);
                
                
                detectWeight_cells = coeffAllCells(:,1:nPCs)*detectWeight;
                targetWeight_cells = coeffAllCells(:,1:nPCs)*targetWeight;

                decodeAnalysis(iexp).av(iav).dvDetect = dv_detect;
                decodeAnalysis(iexp).av(iav).dvTarget = dv_target;
                decodeAnalysis(iexp).av(iav).correlationDetect = detectCorr;
                decodeAnalysis(iexp).av(iav).correlationTarget = targetCorr;
                decodeAnalysis(iexp).av(iav).weightDetect = detectWeight_cells;
                decodeAnalysis(iexp).av(iav).weightTarget = targetWeight_cells;
                decodeAnalysis(iexp).av(iav).weightDetect_pcs = detectWeight;
                decodeAnalysis(iexp).av(iav).weightTarget_pcs = targetWeight;
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
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_distOnly_holdout = pctCorrectDetect_ho_dist;
                decodeAnalysis(iexp).av(iav).weightDetect_distOnly = detectWeight_cells_distOnly;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_tarOnly_holdout = pctCorrectDetect_ho_tar;
                decodeAnalysis(iexp).av(iav).weightDetect_tarOnly = detectWeight_cells_tarOnly;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_shuff_holdout = pctCorrectDetect_ho_shuff;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_shuff_holdout = pctCorrectTarget_ho_shuff;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_withStim_holdout = pctCorrectDetect_ho_withStim;
                decodeAnalysis(iexp).av(iav).weightDetect_withStim = detectWeight_cells_withStim;
                decodeAnalysis(iexp).av(iav).modeledStimWeight = modeledStimWeight;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_withChoice_holdout = pctCorrectTarget_ho_withChoice;
                decodeAnalysis(iexp).av(iav).weightTarget_withChoice = targetWeight_cells_withChoice;
                decodeAnalysis(iexp).av(iav).modeledChoiceWeight = modeledChoiceWeight;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_withFixedStim_holdout = pctCorrectDetect_ho_withFixedStim;
                decodeAnalysis(iexp).av(iav).weights_fixedStimGivenChoice_DT = [min(dv_stimEaTrial) max(dv_stimEaTrial)];
                decodeAnalysis(iexp).av(iav).weightDetect_withFixedStim = coeffAllCells(:,1:nPCs)*detectGLM_withFixedStim.beta(2:end-1);
                decodeAnalysis(iexp).av(iav).modeledStimWeight_fixed = detectGLM_withFixedStim.beta(end);
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_withFixedChoice_holdout = pctCorrectTarget_ho_withFixedChoice;
                decodeAnalysis(iexp).av(iav).weights_fixedChoiceGivenStim_NY = [min(dv_choiceEaTrial) max(dv_choiceEaTrial)];
                decodeAnalysis(iexp).av(iav).weightTarget_withFixedChoice = coeffAllCells(:,1:nPCs)*targetGLM_withFixedChoice.beta(2:end-1);
                decodeAnalysis(iexp).av(iav).modeledChoiceWeight_fixed = targetGLM_withFixedChoice.beta(end);
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_withFixedStim_holdout_pcsOnly = ...
                    pctCorrectDetect_ho_withFixedStim_pcsOnly;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_withFixedChoice_holdout_pcsOnly = ...
                    pctCorrectTarget_ho_withFixedChoice_pcsOnly;
                
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_choiceOnly = pctCorrectTarget_choiceOnly;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_stimOnly = pctCorrectDetect_stimOnly;
                
%                 decodeAnalysis(iexp).av(iav).choiceModel_fixedStim_bootStruct = choiceModel_fixedStim_struct;
%                 decodeAnalysis(iexp).av(iav).stimModel_fixedChoice_bootStruct = stimModel_fixedChoice_struct;
                
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_resid = pctCorrectTarget_ho_resid;
                decodeAnalysis(iexp).av(iav).weightTarget_resid = targetWeight_cells_resid;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_resid = pctCorrectDetect_ho_resid;
                decodeAnalysis(iexp).av(iav).weightDetect_resid = detectWeight_cells_resid;
                
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
%                         catchTargetResp = catchResp(:,cellInd);
                        catchTargetResp = catchResp(:,1:nPCs);
                        nCatch = length(catchDetectInd);

                            [catchPctCorrectDetect,catchCorrectTrialDetect] = getPctCorr_trainData(...
                                detectGLMOther{iav},catchTargetResp,catchDetectInd,...
                                decodeAnalysis(iexp).av(iav).dvDetect);
                            [catchPctCorrectTarget,catchCorrectTrialTarget] = getPctCorr_trainData(...
                                targetGLMOther{iav},catchTargetResp,catchTargetInd,...
                                decodeAnalysis(iexp).av(iav).dvTarget);
                            [catchPctCorrectDetect_aud,catchCorrectTrialDetect_aud] = getPctCorr_trainData(...
                                detectGLMOther{auditoryTrials},catchTargetResp,...
                                catchDetectInd,...
                                decodeAnalysis(iexp).av(auditoryTrials).dvDetect);            
                            [catchPctCorrectTarget_aud,catchCorrectTrialTarget_aud] = getPctCorr_trainData(...
                                targetGLMOther{auditoryTrials},catchTargetResp,...
                                catchTargetInd,...
                                decodeAnalysis(iexp).av(auditoryTrials).dvTarget);
%                             catchPctCorrectDetect = getPctCorr_trainData(...
%                                 detectGLMOther{iav},catchRespBalanced,catchDetectIndBalanced,...
%                                 decodeAnalysis(iexp).av(iav).dvDetect);
%                             catchPctCorrectTarget = getPctCorr_trainData(...
%                                 targetGLMOther{iav},catchRespBalanced,catchTargetIndBalanced,...
%                                 decodeAnalysis(iexp).av(iav).dvTarget);
%                             catchPctCorrectDetect_audModel = getPctCorr_trainData(...
%                                 detectGLMOther{auditoryTrials},catchRespBalanced,...
%                                 catchDetectIndBalanced,...
%                                 decodeAnalysis(iexp).av(auditoryTrials).dvDetect);            
%                             catchPctCorrectTarget_audModel = getPctCorr_trainData(...
%                                 targetGLMOther{auditoryTrials},catchRespBalanced,...
%                                 catchTargetIndBalanced,...
%                                 decodeAnalysis(iexp).av(auditoryTrials).dvTarget);
%                         else
%                             catchPctCorrectDetect = nan;
%                             catchCorrectTrialDetect = nan(0,1);
%                             catchPctCorrectTarget = nan;
%                             catchCorrectTrialTarget = nan(0,1);
%                             catchPctCorrectDetect_audModel = nan;
%                             catchCorrectTrialDetect_audModel = nan(0,1);        
%                             catchPctCorrectTarget_audModel = nan;
%                             catchCorrectTrialTarget_audModel = nan(0,1);  
%                         end
                        

                        trStimID = discretize(decodeDataExpt(iexp).av(iav).stim,oriBins);
                        matchedTrialsID = trStimID(trSampleIndOther{iav});
                        resp_val = respOther{iav};
%                         catchStimID = cat(2,ones(1,nCatch),...
%                             discretize(decodeDataExpt(iexp).av(iav).catchStim,oriBins));
                        catchStimID = discretize(decodeDataExpt(iexp).av(iav).catchStim,oriBins);
                        
                        [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{iav});
                        if ~any(ismember(matchedTrialsID,unique(catchStimID)))
                            catchMatchTrials = [];
                        else
                            catchMatchTrials = cell2mat(getMatchedValidTrialIndex(...
                                matchedTrialsID,catchStimID)); 
                        end
                        [validCatchMatchDetect,validCatchCorrectTrialDetect] = getPctCorr_hoData_subGroup(...
                            resp_val,detectTrInd,catchMatchTrials,decodeAnalysis(iexp).av(iav).dvDetect);  
                        [validCatchMatchTarget,validCatchCorrectTrialTarget] = getPctCorr_hoData_subGroup(...
                            resp_val,targetTrInd,catchMatchTrials,decodeAnalysis(iexp).av(iav).dvTarget);
                        
%                         trStimID = discretize(decodeDataExpt(iexp).av(auditoryTrials).stim,ampBins);
%                         matchedTrialsID = trStimID(trSampleIndOther{auditoryTrials});
%                         matchedTrialsID(matchedTrialsID > 1) = 2;
%                         catchStimID(catchStimID > 1) = 2;
%                         [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{auditoryTrials});
%                         catchMatchTrials = cell2mat(getMatchedValidTrialIndex(...
%                             matchedTrialsID,catchStimID));
                        [validCatchMatchDetect_aud,validCatchCorrectTrialDetect_aud] = getPctCorr_trainData(...
                            detectGLMOther{otherAV},resp_val(catchMatchTrials,:),...
                            detectTrInd(catchMatchTrials),decodeAnalysis(iexp).av(otherAV).dvDetect);
                        [validCatchMatchTarget_aud,validCatchCorrectTrialTarget_aud] = getPctCorr_trainData(...
                            targetGLMOther{otherAV},resp_val(catchMatchTrials,:),...
                            targetTrInd(catchMatchTrials),decodeAnalysis(iexp).av(otherAV).dvTarget); 
                       
                    else
                        validCatchMatchDetect = [];
                        validCatchCorrectTrialDetect = [];
                        validCatchMatchTarget = [];
                        validCatchCorrectTrialTarget = [];
                        validCatchMatchDetect_aud = [];
                        validCatchCorrectTrialDetect_aud = [];
                        validCatchMatchTarget_aud = [];
                        validCatchCorrectTrialTarget_aud = [];
                        catchPctCorrectDetect = [];
                        catchPctCorrectTarget = [];
                        catchPctCorrectDetect_aud = [];
                        catchPctCorrectTarget_aud = [];
                        catchCorrectTrialDetect = [];
                        catchCorrectTrialTarget = [];
                        catchCorrectTrialDetect_aud = [];
                        catchCorrectTrialTarget_aud = [];
                        
                    end
                elseif iav == 2
                    otherAV = 1;
                    validCatchMatchDetect = [];
                    validCatchMatchTarget = [];
                    validCatchMatchDetect_aud = [];
                    validCatchMatchTarget_aud = [];
                    catchPctCorrectDetect = [];
                    catchPctCorrectTarget = [];
                    catchPctCorrectDetect_aud = [];
                    catchPctCorrectTarget_aud = [];
                end
                resp = respOther{otherAV};
                [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOutOther{otherAV});

%                 dv_detect = mean(detectTrInd);
%                 dv_target = mean(targetTrInd);
                [detecttrind4model,targettrind4model] = getStimAndBehaviorYs(trOutOther{iav});
                dv_detect = mean(detecttrind4model);
                dv_target = mean(targettrind4model);
                
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
                
                % test other model with fixedstim/choice predictors
                pctCorrectDetect_withFixedStim = getPctCorr_weightsOnly(...
                    detectGLMOther_withStim{iav}.beta(1:(nPCs+1)),resp,...
                    detectTrInd,dv_detect);
                pctCorrectTarget_withFixedChoice = getPctCorr_weightsOnly(...
                    targetGLMOther_withChoice{iav}.beta(1:(nPCs+1)),resp,...
                    targetTrInd,dv_target);
                
                % test other distractor only detect model
                resp = respOther_dist{otherAV};
                [detectTrInd, ~] = getStimAndBehaviorYs(trOutOther_dist{otherAV});
                dv_detect = mean(getStimAndBehaviorYs(trOutOther_dist{iav}));
                pctCorrectDetect_dist = getPctCorr_trainData(detectGLMOther_dist{iav},resp,detectTrInd,dv_detect);
                % test other target only detect model
                resp = respOther_tar{otherAV};
                [detectTrInd, ~] = getStimAndBehaviorYs(trOutOther_tar{otherAV});
                dv_detect = mean(getStimAndBehaviorYs(trOutOther_tar{iav}));
                pctCorrectDetect_tar = getPctCorr_trainData(detectGLMOther_tar{iav},resp,detectTrInd,dv_detect);
                
                % test other shuffled model
                [resp_detect,resp_target] = deal(respOther_shuff{otherAV}{:});
                [detectTrInd,targetTrInd] = deal(trOutInd_shuff{otherAV}{:});
                dv_detect = mean(trOutInd_shuff{iav}{1});
                dv_target = mean(trOutInd_shuff{iav}{2});    
                pctCorrectDetect_shuff = getPctCorr_trainData(...
                    detectGLMOther_shuff{iav},resp_detect,detectTrInd,dv_detect);
                pctCorrectTarget_shuff = getPctCorr_trainData(...
                    targetGLMOther_shuff{iav},resp_target,targetTrInd,dv_target);
                
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_otherAV = pctCorrectDetect;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_otherAV = pctCorrectTarget;
                decodeAnalysis(iexp).av(iav).pctCorrectDetectxStim_otherAV = pctCorrectDetectxStim;
                decodeAnalysis(iexp).av(iav).pctCorrectTargetxStim_otherAV = pctCorrectTargetxStim;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectDetect_holdout = validCatchMatchDetect;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectTarget_holdout = validCatchMatchTarget;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectDetect = catchPctCorrectDetect;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectTarget = catchPctCorrectTarget;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectDetect_testAudModel = catchPctCorrectDetect_aud;
                decodeAnalysis(iexp).av(iav).invalidPctCorrectTarget_testAudModel = catchPctCorrectTarget_aud;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectDetect_testAudModel = validCatchMatchDetect_aud;
                decodeAnalysis(iexp).av(iav).validMatchedPctCorrectTarget_testAudModel = validCatchMatchTarget_aud;
                
                decodeAnalysis(iexp).av(iav).validMatchedCorrectTrialsDetect = validCatchCorrectTrialDetect;
                decodeAnalysis(iexp).av(iav).validMatchedCorrectTrialsTarget = validCatchCorrectTrialTarget;
                decodeAnalysis(iexp).av(iav).invalidCorrectTrialsDetect = catchCorrectTrialDetect;
                decodeAnalysis(iexp).av(iav).invalidCorrectTrialsTarget = catchCorrectTrialTarget;
                decodeAnalysis(iexp).av(iav).invalidCorrectTrialsDetect_testAud = catchCorrectTrialDetect_aud;
                decodeAnalysis(iexp).av(iav).invalidCorrectTrialsTarget_testAud = catchCorrectTrialTarget_aud;
                decodeAnalysis(iexp).av(iav).validMatchedCorrectTrialsDetect_testAud = validCatchCorrectTrialDetect_aud;
                decodeAnalysis(iexp).av(iav).validMatchedCorrectTrialsTarget_testAud = validCatchCorrectTrialTarget_aud;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_distOnly_otherAV = pctCorrectDetect_dist;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_tarOnly_otherAV = pctCorrectDetect_tar;
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_shuff_otherAV = pctCorrectDetect_shuff;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_shuff_otherAV = pctCorrectTarget_shuff;
                
                decodeAnalysis(iexp).av(iav).pctCorrectDetect_withFixedStim_otherAV = ...
                    pctCorrectDetect_withFixedStim;
                decodeAnalysis(iexp).av(iav).pctCorrectTarget_withFixedChoice_otherAV = ...
                    pctCorrectTarget_withFixedChoice;
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
            
            
            decodeAnalysis(iexp).comboTrainWeightDetect = coeffAllCells(:,1:nPCs)*detectGLM.beta(2:end);
            decodeAnalysis(iexp).comboTrainWeightTarget = coeffAllCells(:,1:nPCs)*targetGLM.beta(2:end);
            decodeAnalysis(iexp).av(visualTrials).pctCorrectDetect_comboTrain = pctCorrectDetect_combo_vis;
            decodeAnalysis(iexp).av(visualTrials).pctCorrectTarget_comboTrain = pctCorrectTarget_combo_vis;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectDetect_comboTrain = pctCorrectDetect_combo_aud;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectTarget_comboTrain = pctCorrectTarget_combo_aud;
            decodeAnalysis(iexp).av(visualTrials).pctCorrectDetect_visTrainComboMatch = pctCorrectDetect_vis;
            decodeAnalysis(iexp).av(visualTrials).pctCorrectTarget_visTrainComboMatch = pctCorrectTarget_vis;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectDetect_audTrainComboMatch = pctCorrectDetect_aud;
            decodeAnalysis(iexp).av(auditoryTrials).pctCorrectTarget_audTrainComboMatch = pctCorrectTarget_aud;
        end
end
    
%% plot the performance of the choice model on all data vs. data - stimulus response

dc_attn = decodeAnalysis;
avName = {'Vis','Aud'};
modelName = {'Stimulus','Choice'};
nexp_attn = size(dc_attn,2);
stimModInd = 1;
choiceModInd = 2;

%model performance - attention mice
mdl_attn = struct;
for iav = 1:2
    mdl_attn(iav).name = avName{iav};
    mdl_attn(iav).mdl(stimModInd).name = 'Stimulus';
    mdl_attn(iav).mdl(choiceModInd).name = 'Choice';
    mdl_attn(iav).mdl(stimModInd).pctCorr_all = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_all = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_all(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectAllTarget_holdout;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_all(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectAllDetect_holdout;
    end
end
goodPerfInd_attn = mdl_attn(visualTrials).mdl(stimModInd).pctCorr_all...
    > pctCorrThresh;
% choice model with stimulus responses subtracted
for iav = 1:2
    mdl_attn(iav).mdl(stimModInd).pctCorr_all_resid = nan(1,nexp_attn);
    mdl_attn(iav).mdl(choiceModInd).pctCorr_all_resid = nan(1,nexp_attn);
    for iexp = 1:nexp_attn
        mdl_attn(iav).mdl(stimModInd).pctCorr_all_resid(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectTarget_resid;
        mdl_attn(iav).mdl(choiceModInd).pctCorr_all_resid(iexp) = ...
            dc_attn(iexp).av(iav).pctCorrectDetect_resid;
    end
    
end

figure
suptitle('Expts with good perf on vis stim model')
for iav = 1:2
    for imod = 1:2
        if imod == 1
            iplot = 0;
        else
            iplot = 2;
        end
        subplot(2,2,iav+iplot)
        y = [mdl_attn(iav).mdl(imod).pctCorr_all(goodPerfInd_attn)',...
            mdl_attn(iav).mdl(imod).pctCorr_all_resid(goodPerfInd_attn)'];
        plot(1:2,y,'k-')
        hold on
        errorbar(1:2,mean(y,1),ste(y,1),'.')
        [~,p] = ttest(y(:,1),y(:,2));
        title(sprintf('Attn %s %s Model, p=%s',mdl_attn(iav).name,...
            mdl_attn(iav).mdl(imod).name, sigfigString(p)))
        figXAxis([],'',[0 3],1:2,{'All Data', 'Data - (mean stim resp)'})
        ax = gca;
        ax.XTickLabelRotation = -45;
        figYAxis([],'Pct Corr.',[0 1])
        figAxForm
        hline(pctCorrThresh,'k:')
    end
end