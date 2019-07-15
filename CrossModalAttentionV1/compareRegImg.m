clear all
close all
ds = 'FSAV_V1_SOM_naive';
rc = behavConstsAV;
eval(ds)
slct_expt = 2;
%%
for iexp = slct_expt

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;

%% load all behavior data
% concatenate data files
for irun = 1:expt(iexp).nrun

    if isfield(expt,'nframesPerRun')
        if ~isempty(expt(iexp).nframesPerRun)
            nFr = expt(iexp).nframesPerRun(irun);
        else
            nFr = [];
        end
    else
        nFr = [];
    end
    
    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    if isfield (expt,'greenredsimultaneous')
        if expt(iexp).greenredsimultaneous == 1
            input = loadMworksFile(SubNum,expDate,expTime);
            data_temp = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName);
        else
            if ~isempty(nFr)
                [input, data_temp, t] = Load_SBXdataPlusMWorksData(...
                    SubNum,expDate,expTime,mouse,runFolder,fName,nFr);
                disp(t)
            else
                [input, data_temp, t] = Load_SBXdataPlusMWorksData(...
                    SubNum,expDate,expTime,mouse,runFolder,fName);
                disp(t)
            end
        end
    else
        if ~isempty(nFr)
            [input, data_temp, t] = Load_SBXdataPlusMWorksData(...
                SubNum,expDate,expTime,mouse,runFolder,fName);
            disp(t)
        else
            [input, data_temp, t] = Load_SBXdataPlusMWorksData(...
                SubNum,expDate,expTime,mouse,runFolder,fName,nFr);
            disp(t)
        end
        
    end
    
    if irun == 1
        data_bx = data_temp;
        input_bx = input;
    else
        data_bx = cat(3,data_bx,data_temp);
        try
            input_bx = [input_bx input];
        catch
            inpNames1 = fieldnames(input_bx);
            inpNames2 = fieldnames(input);
            inpLong = gt(length(inpNames1),length(inpNames2));
            if inpLong == 1
                inpPlusInd = ismember(inpNames1,inpNames2);
                inpPlus = inpNames1(~inpPlusInd);
                for i = 1:length(inpPlus)
                    input.(genvarname(inpPlus{i})) = cell(1,input.trialSinceReset);
                end
            else
                inpPlusInd = ismember(inpNames2,inpNames1);
                inpPlus = inpNames2(~inpPlusInd);
                for i = 1:length(inpPlus)
                    input_bx.(char(genvarname(inpPlus(i)))) = cell(1,length(inpNames1));
                end
            end
            input_temp = [input_bx input];
        end
    end
    clear data_temp input
end
input_bx = concatenateDataBlocks(input_bx);
nfr_bx = size(data_bx,3);

% remove negative data by subtraction
data_sub = data_bx-min(min(min(data_bx,[],1),[],2),[],3);
data_bx = data_sub;
clear data_sub

%%
fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);

nMeanImages = 8;
[nRows,nCols] = optimizeSubplotDim(nMeanImages+1);

randStarterFrames = sort(randsample(nfr_bx-100,nMeanImages));
setFigParams4Print('landscape')
figure
suptitle([SubNum '-' expDate])
for iimg = 1:nMeanImages
    subplot(nRows,nCols,iimg)
    imagesc(mean(data_bx(:,:,randStarterFrames(iimg):randStarterFrames(iimg)+100),3))
    title(sprintf('%s:%s',num2str(randStarterFrames(iimg)),num2str(randStarterFrames(iimg)+100)))
end
% if doCorrImg
%     try
%         load(fullfile(fnout,'regOuts&Img.mat'));
%         subplot(nRows,nCols,nMeanImages+1)
%         imagesc(data_corr_img)
%         title('corr image')
%     catch
%         disp('no corr image')
%     end
% end
if ~exist(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds),'dir')
    mkdir(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds))
end
print(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds,...
    ['rand image samples_' SubNum '-' expDate]),'-dpdf','-fillpage')
savefig(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds,...
    ['rand image samples_' SubNum '-' expDate]))
clear data_bx
end
clear all