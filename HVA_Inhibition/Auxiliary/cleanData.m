function new_data = cleanData(dStimTrials,data_dim,stimOn,varargin)
p = inputParser;
p.addParamValue('pclamp', true, @islogical);
p.addParamValue('slice', false, @islogical);
p.addParamValue('rebaseline', false, @islogical);
p.addParamValue('smooth',false,@islogical);
p.addParamValue('baselineS', 0.1, @isnumeric);
p.addParamValue('afterStimS', 0.3, @isnumeric);
p.addParamValue('sampleFreq', 20000, @isnumeric);


% parse inputs
p.parse(varargin{:});
params = p.Results;

if params.pclamp
    
    if numel(stimOn) == 1
        nTrials = size(dStimTrials,3);
        stimOn = repmat(stimOn,1,nTrials);
    else
        nTrials = numel(stimOn);
    end
    for trial = 1:nTrials
        new_data.baseline_pts(:,trial) = (squeeze(dStimTrials((stimOn(trial)-params.baselineS*params.sampleFreq):stimOn(trial),data_dim,trial)));
        new_data.baseline(trial) = nanmean(new_data.baseline_pts(:,trial)); % get baseline for 1st pulse


        
        if params.smooth
            new_data.pts(:,trial) = smooth(squeeze(dStimTrials(stimOn(trial):stimOn(trial)+params.afterStimS*params.sampleFreq,data_dim,trial))-new_data.baseline(trial),200);   
        else
            new_data.pts(:,trial) = squeeze(dStimTrials(stimOn(trial):stimOn(trial)+params.afterStimS*params.sampleFreq,data_dim,trial))-new_data.baseline(trial);   
        end
    %     if needs to be rebaselined (large fluctuation that is skewed when
    %     taking mean)
    %     

    if params.rebaseline
        if params.slice
            if mean(new_data.pts(:,trial))>0
                if abs(mean(new_data.pts(50:65,trial))) > 5
                    new_data.pts(:,trial) = new_data.pts(:,trial) - mean(new_data.pts(50:65,trial));
                end 
            else
                if abs(mean(new_data.pts(35:50,trial))) > 5
                    new_data.pts(:,trial) = new_data.pts(:,trial) - mean(new_data.pts(35:50,trial));
                end 
            end
        else
            if abs(mean(new_data.pts(750,trial))) > 20
                new_data.pts(:,trial) = new_data.pts(:,trial) - mean(new_data.pts(500:1000,trial));
            end
        end

    end
    end
else
    new_data.pts = cell2mat(arrayfun(@(x) dStimTrials(x:floor(x+params.afterStimS*params.sampleFreq))-nanmean(dStimTrials(floor(x-params.baselineS*params.sampleFreq):x)),stimOn,'un',0));
    if params.rebaseline
        if any(mean(new_data.pts(750,:)) > 50)
            replace_i = mean(new_data.pts(750,:)) > 50;
            new_data.pts(:,replace_i) = new_data.pts(:,mean(new_data.pts(750,:)) > 50) - nanmean(new_data.pts(250:750,replace_i));
        end

    end
end
