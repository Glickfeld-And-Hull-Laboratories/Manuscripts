% attend = 1;
cc(1,:)=[0 0 0];
cc(2,:)=[0.75 0 0.75];
attn_label = {'attn';'no attn'};
av_label = {'vis';'aud'};
number=false;
ssize=20;

mi_lim = [0 1];

%% get PID info from results structure

if doBootstrapPID
    PID = struct;
    PID.II.Vis = nan(nBoot,length(D));
    PID.MI_SR.Vis = nan(nBoot,length(D));
    PID.MI_BR.Vis = nan(nBoot,length(D));
    PID.MI_SB.Vis = nan(nBoot,length(D));
    PID.II.Aud = nan(nBoot,length(D));
    PID.MI_SR.Aud = nan(nBoot,length(D));
    PID.MI_BR.Aud = nan(nBoot,length(D));
    PID.MI_SB.Aud = nan(nBoot,length(D));
    for iboot = 1:nBoot
        PID.II.Vis(iboot,:) = PIDresults.PID{iboot}.II.Vis;
        PID.MI_SR.Vis(iboot,:) = PIDresults.PID{iboot}.MI_SR.Vis;
        PID.MI_BR.Vis(iboot,:) = PIDresults.PID{iboot}.MI_BR.Vis;
        PID.MI_SB.Vis(iboot,:) = PIDresults.PID{iboot}.MI_SB.Vis;
        PID.II.Aud(iboot,:) = PIDresults.PID{iboot}.II.Aud;
        PID.MI_SR.Aud(iboot,:) = PIDresults.PID{iboot}.MI_SR.Aud;
        PID.MI_BR.Aud(iboot,:) = PIDresults.PID{iboot}.MI_BR.Aud;
        PID.MI_SB.Aud(iboot,:) = PIDresults.PID{iboot}.MI_SB.Aud;      
    end
else
    PID = PIDresults.PID;
end
%% plot PID results across experiments
if ~doMultiPC
setFigParams4Print('landscape')
fig_pid = figure;
fig_mi = figure;
if doBootstrapPID
for iattend = 0:1
    if iattend == 0
        h = [];
    end
    idx = cell2mat({data.hasAttention}) == iattend;
    figure(fig_mi)
    subplot 131; hold on
    x = mean(PID.MI_SR.Vis(:,idx),1);
    [xerr_l,xerr_u] = ciFromBoot(PID.MI_SR.Vis(:,idx),95);
    y = mean(PID.MI_SR.Aud(:,idx),1);
    [yerr_l,yerr_u] = ciFromBoot(PID.MI_SR.Aud(:,idx),95);
    errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'.','MarkerSize',ssize)
%     plot(x,y,'.','MarkerSize',ssize)
    title('mutual info in the stimulus and response')
    subplot 132; hold on
    x = mean(PID.MI_BR.Vis(:,idx),1);
    [xerr_l,xerr_u] = ciFromBoot(PID.MI_BR.Vis(:,idx),95);
    y = mean(PID.MI_BR.Aud(:,idx),1);
    [yerr_l,yerr_u] = ciFromBoot(PID.MI_BR.Aud(:,idx),95);
    h(iattend+1) = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'.','MarkerSize',ssize);
    if iattend == 1
        legend(attn_label,'location','northeast')
        refline(1)
        for iplot = 1:3
        subplot(1,3,iplot)
        figXAxis([],'Visual',mi_lim)
        figYAxis([],'Auditory',mi_lim)
        figAxForm
        refline(1)
        end
    end
%     plot(x,y,'.','MarkerSize',ssize)
%     plot(PID.MI_BR.Vis(idx),PID.MI_BR.Aud(idx),'.','MarkerSize',ssize)
    title('mutual info in the choice and response')
    subplot 133; hold on
    x = mean(PID.MI_SB.Vis(:,idx),1);
    [xerr_l,xerr_u] = ciFromBoot(PID.MI_SB.Vis(:,idx),95);
    y = mean(PID.MI_SB.Aud(:,idx),1);
    [yerr_l,yerr_u] = ciFromBoot(PID.MI_SB.Aud(:,idx),95);
    errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'.','MarkerSize',ssize);
%     plot(x,y,'.','MarkerSize',ssize)
%     h(iattend+1) = plot(PID.MI_SB.Vis(idx),PID.MI_SB.Aud(idx),'.','MarkerSize',ssize);
    title('mutual info in the stimulus and choice')
    
    
    figure(fig_pid)
    subplot 321; hold on
    x = mean(PID.II.Vis(:,idx)./PID.MI_SR.Vis(:,idx),1);
    [xerr_l,xerr_u] = ciFromBoot(PID.II.Vis(:,idx)./PID.MI_SR.Vis(:,idx),95);
    y = mean(PID.II.Aud(:,idx)./PID.MI_SR.Aud(:,idx),1);
    [yerr_l,yerr_u] = ciFromBoot(PID.II.Aud(:,idx)./PID.MI_SR.Aud(:,idx),95);
    errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'.','MarkerSize',ssize)
%     plot(x,y,'.','MarkerSize',ssize)
%     scatter(PID.II.Vis(idx)./PID.MI_SR.Vis(idx),PID.II.Aud(idx)./PID.MI_SR.Aud(idx),ssize)
    axis([0 1 0 1])
    title('Percent Sensory info in Pop that drives behavior  (PID/SR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm
    if iattend == 1
        legend(attn_label,'location','southwestoutside')
        refline(1)
    end

    subplot 322; hold on
    x = mean(1-PID.II.Vis(:,idx)./PID.MI_SR.Vis(:,idx),1);
    [xerr_l,xerr_u] = ciFromBoot(1-PID.II.Vis(:,idx)./PID.MI_SR.Vis(:,idx),95);
    y = mean(1-PID.II.Aud(:,idx)./PID.MI_SR.Aud(:,idx),1);
    [yerr_l,yerr_u] = ciFromBoot(1-PID.II.Aud(:,idx)./PID.MI_SR.Aud(:,idx),95);
    errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'.','MarkerSize',ssize)
%     plot(x,y,'.','MarkerSize',ssize)    
%     scatter(1-PID.II.Vis(idx)./PID.MI_SR.Vis(idx),1-PID.II.Aud(idx)./PID.MI_SR.Aud(idx),ssize)
    axis([0 1 0 1])
    refline(1)
    title('Percent Sensory info in Pop that does not drives behavior (1-PID/SR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm

    subplot 323; hold on
    x = mean(1-PID.II.Vis(:,idx)./PID.MI_BR.Vis(:,idx),1);
    [xerr_l,xerr_u] = ciFromBoot(1-PID.II.Vis(:,idx)./PID.MI_BR.Vis(:,idx),95);
    y = mean(1-PID.II.Aud(:,idx)./PID.MI_BR.Aud(:,idx),1);
    [yerr_l,yerr_u] = ciFromBoot(1-PID.II.Aud(:,idx)./PID.MI_BR.Aud(:,idx),95);
    errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'.','MarkerSize',ssize)
%     plot(x,y,'.','MarkerSize',ssize)    
%     scatter(1-PID.II.Vis(idx)./PID.MI_BR.Vis(idx),1-PID.II.Aud(idx)./PID.MI_BR.Aud(idx),ssize)
    axis([0 1 0 1])
    refline(1)
    title('Percent Choice info in Pop not related to Stim (1-PID/BR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm

    subplot 324; hold on
    x = mean(1-PID.II.Vis(:,idx)./PID.MI_SB.Vis(:,idx),1);
    [xerr_l,xerr_u] = ciFromBoot(1-PID.II.Vis(:,idx)./PID.MI_SB.Vis(:,idx),95);
    y = mean(1-PID.II.Aud(:,idx)./PID.MI_SB.Aud(:,idx),1);
    [yerr_l,yerr_u] = ciFromBoot(1-PID.II.Aud(:,idx)./PID.MI_SB.Aud(:,idx),95);
    errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'.','MarkerSize',ssize)
%     plot(x,y,'.','MarkerSize',ssize)    
%     scatter(1-PID.II.Vis(idx)./PID.MI_SB.Vis(idx),1-PID.II.Aud(idx)./PID.MI_SB.Aud(idx),ssize)
    axis([0 1 0 1])
    refline(1)
    title('Percent Stimulus info that drives behavior via an alt. pathway (1-PID/SB)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm

%     figure(fig_m_pid)
    subplot(3,2,iattend+5);hold on

    % non_readout_sensory_info=S_R_info-I_II;
    d = 1-[mean(PID.II.Vis(:,idx)./PID.MI_SR.Vis(:,idx),1)',mean(PID.II.Aud(:,idx)./PID.MI_SR.Aud(:,idx),1)'];
    b1=mean(d); % sensory information in neural response R that is not read out for behavior
    b1err = ste(d,1);
    

    %internal_choice_info=C_R_info-I_II;
    d = 1-[mean(PID.II.Vis(:,idx)./PID.MI_BR.Vis(:,idx))',mean(PID.II.Aud(:,idx)./PID.MI_BR.Aud(:,idx))'];
    b2=mean(d);% choice information in neural response R that is not related to the stimulus
    b2err = ste(d,1);

    %S_C_info_from_unobserved_R=S_C_info-I_II;
    d = 1-[mean(PID.II.Vis(:,idx)./PID.MI_SB.Vis(:,idx))',mean(PID.II.Aud(:,idx)./PID.MI_SB.Aud(:,idx))'];
    b3=mean(d);% the part of "behavioral performance" that cannot be explained with recorded neural feature R
    b3err = ste(d,1);
        
%     bar([b1;b2;b3;]) 
    for iav = 1:2
        h = errorbar([b1(iav);b2(iav);b3(iav);],[b1err(iav);b2err(iav);b3err(iav);],'.','MarkerSize',20);
        h.Color = cc(iav,:);
    end
    
    legend(av_label,'location','northeastoutside')
    figXAxis([],'',[0 4],1:3,{'1-II/SR';'1-II/BR';'1-II/SB'})
    xtickangle(-45)
    figYAxis([],'II',[0 1])
    figAxForm
    title({attn_label{iattend+1};'error across expt'})
end
else
for iattend = 0:1
    if iattend == 0
        h = [];
    end
    idx = cell2mat({data.hasAttention}) == iattend;
    figure(fig_pid)
    subplot 431; hold on
    plot(PID.MI_SR.Vis(idx),PID.MI_SR.Aud(idx),'.','MarkerSize',ssize)
    title('mutual info in the stimulus and response')
    subplot 432; hold on
    plot(PID.MI_BR.Vis(idx),PID.MI_BR.Aud(idx),'.','MarkerSize',ssize)
    title('mutual info in the choice and response')
    subplot 433; hold on
    h(iattend+1) = plot(PID.MI_BR.Vis(idx),PID.MI_BR.Aud(idx),'.','MarkerSize',ssize);
    title('mutual info in the stimulus and choice')
    if iattend == 1
        legend(h,attn_label)
    end
    for iplot = 1:3
        subplot(4,3,iplot)
        figXAxis([],'Visual',mi_lim)
        figYAxis([],'Auditory',mi_lim)
        figAxForm
        refline(1)
    end
    
%     figure(fig_pid)
    subplot 423; hold on
    scatter(PID.II.Vis(idx)./PID.MI_SR.Vis(idx),PID.II.Aud(idx)./PID.MI_SR.Aud(idx),ssize)
    axis([0 1 0 1])
    title('Percent Sensory info in Pop that drives behavior  (PID/SR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm
    if iattend == 1
        legend(attn_label,'location','southwestoutside')
        refline(1)
    end

    subplot 424; hold on
    scatter(1-PID.II.Vis(idx)./PID.MI_SR.Vis(idx),1-PID.II.Aud(idx)./PID.MI_SR.Aud(idx),ssize)
    axis([0 1 0 1])
    refline(1)
    title('Percent Sensory info in Pop that does not drives behavior (1-PID/SR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm

    subplot 425; hold on
    scatter(1-PID.II.Vis(idx)./PID.MI_BR.Vis(idx),1-PID.II.Aud(idx)./PID.MI_BR.Aud(idx),ssize)
    axis([0 1 0 1])
    refline(1)
    title('Percent Choice info in Pop not related to Stim (1-PID/BR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm

    subplot 426; hold on
    scatter(1-PID.II.Vis(idx)./PID.MI_SB.Vis(idx),1-PID.II.Aud(idx)./PID.MI_SB.Aud(idx),ssize)
    axis([0 1 0 1])
    refline(1)
    title('Percent Stimulus info that drives behavior via an alt. pathway (1-PID/SB)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm

%     figure(fig_m_pid)
    subplot(4,2,iattend+7)

    % non_readout_sensory_info=S_R_info-I_II;
    b1=mean([(PID.II.Vis(idx)./PID.MI_SR.Vis(idx))',(PID.II.Aud(idx)./PID.MI_SR.Aud(idx))']);
    b1=1-b1;% sensory information in neural response R that is not read out for behavior

    %internal_choice_info=C_R_info-I_II;
    b2=mean([(PID.II.Vis(idx)./PID.MI_BR.Vis(idx))',(PID.II.Aud(idx)./PID.MI_BR.Aud(idx))']);
    b2=1-b2;% choice information in neural response R that is not related to the stimulus

    %S_C_info_from_unobserved_R=S_C_info-I_II;
    b3=mean([(PID.II.Vis(idx)./PID.MI_SB.Vis(idx))',(PID.II.Aud(idx)./PID.MI_SB.Aud(idx))']);
    b3 = 1-b3; % the part of "behavioral performance" that cannot be explained 
               % with recorded neural feature R
        
    bar([b1;b2;b3;]) 
    
    figXAxis([],'',[0 4],1:3,{'1-II/SR';'1-II/BR';'1-II/SB'})
    figYAxis([],'II',[0 1])
    figAxForm
    title(attn_label{iattend+1})
    legend(av_label,'location','northwest')
end
end
if doLoadPrevious
    saveStr = ['_results_' runTag '_' prevTimestampID];
else
    saveStr = ['_results_' runTag '_' timestampID];
end
    figure(fig_pid)
    print(fullfile(fnout,['PID' saveStr]),'-dpdf','-fillpage')
    figure(fig_mi)
    print(fullfile(fnout,['MI' saveStr]),'-dpdf','-fillpage')
end

%% compare npcs used
doExptPlot = false;
exptNPCsID = '200721_1946';
mi_lim_R = [0 0.5];
mi_lim_SB = [0 0.75];
if doMultiPC
    load(fullfile(fnout,['PID_results_' exptNPCsID]));
    y_SR = nan(length(data),length(npcs)+1,2);
    y_BR = nan(length(data),length(npcs)+1,2);
    y_SB = nan(length(data),length(npcs)+1,2);
    x = nan(length(data),length(npcs)+1);
    for iexp = 1:length(data)
        if doExptPlot
            figure
            suptitle(data(iexp).expt_name)
        end
        for iav = 0:1
            if iav == 0
                ntrials = length(PIDresults.D{iexp}.Vis.Rlabels);
            else
                ntrials = length(PIDresults.D{iexp}.Aud.Rlabels);
            end
            for ipc = 1:length(npcs)
                if iav == 0
                    y_SR(iexp,ipc,iav+1) = PID{ipc}.MI_SR.Vis(iexp);
                    y_BR(iexp,ipc,iav+1) = PID{ipc}.MI_BR.Vis(iexp);
                    y_SB(iexp,ipc,iav+1) = PID{ipc}.MI_SB.Vis(iexp);
                    if ipc == length(npcs)
                        y_SR(iexp,ipc+1,iav+1) = PIDresults.PID.MI_SR.Vis(iexp);
                        y_BR(iexp,ipc+1,iav+1) = PIDresults.PID.MI_BR.Vis(iexp);
                        y_SB(iexp,ipc+1,iav+1) = PIDresults.PID.MI_SB.Vis(iexp);
                    end
                else
                    y_SR(iexp,ipc,iav+1) = PID{ipc}.MI_SR.Aud(iexp);
                    y_BR(iexp,ipc,iav+1) = PID{ipc}.MI_BR.Aud(iexp);
                    y_SB(iexp,ipc,iav+1) = PID{ipc}.MI_SB.Aud(iexp);
                    if ipc == length(npcs)
                        y_SR(iexp,ipc+1,iav+1) = PIDresults.PID.MI_SR.Aud(iexp);
                        y_BR(iexp,ipc+1,iav+1) = PIDresults.PID.MI_BR.Aud(iexp);
                        y_SB(iexp,ipc+1,iav+1) = PIDresults.PID.MI_SB.Aud(iexp);
                    end
                end
                if ipc == length(npcs)
                    x(iexp,:) = cat(2,npcs,size(PIDresults.D{iexp}.Vis.R,2));
                end
            end
            if doExptPlot
                subplot(2,3,1+(iav*3))
                plot(y_SR(iexp,:,iav+1),'.-','MarkerSize',20,'LineWidth',1)
                title(sprintf('%s Stimulus-Response, %s trials',av_label{iav+1},num2str(ntrials)))
                figXAxis([],'n PCs',[0 length(npcs)+2],1:length(npcs)+1,[npcs,0])
                figYAxis([],'MI (bits)',[])
                figAxForm
                subplot(2,3,2+(iav*3))
                plot(y_BR(iexp,:,iav+1),'.-','MarkerSize',20,'LineWidth',1)
                title(sprintf('%s Behavior-Response, %s trials',av_label{iav+1},num2str(ntrials)))
                figXAxis([],'n PCs',[0 length(npcs)+2],1:length(npcs)+1,[npcs,0])
                figYAxis([],'MI (bits)',[])
                figAxForm
                subplot(2,3,3+(iav*3))
                plot(y_SB(iexp,:,iav+1),'.-','MarkerSize',20,'LineWidth',1)
                title(sprintf('%s Stimulus-Behavior, %s trials',av_label{iav+1},num2str(ntrials)))
                figXAxis([],'n PCs',[0 length(npcs)+2],1:length(npcs)+1,[npcs,0])
                figYAxis([],'MI (bits)',[])
                figAxForm
            end
        end
    end
    figure
    suptitle('nPCs=0 is nPCs=(0.1*ntrials)')
    for iav = 0:1
        subplot(2,3,1+(iav*3)); hold on
        d  = y_SR(:,:,iav+1);
        for iexp = 1:length(data)
            [x_sort,sort_ind] = sort(x(iexp,:));
            plot(x_sort,d(iexp,sort_ind),'k-')
        end
%         plot(d','k-')
        errorbar(mean(sort(x,2)),mean(d,1),ste(d,1),'.-','MarkerSize',20,'LineWidth',1)
        title(sprintf('%s Stimulus-Response',av_label{iav+1}))
        figXAxis([],'n PCs',[0 max(x(:))+1],0:5:max(x(:)),0:5:max(x(:)))
        figYAxis([],'MI (bits)',[mi_lim_R])
        figAxForm
        subplot(2,3,2+(iav*3)); hold on
        d  = y_BR(:,:,iav+1);
        for iexp = 1:length(data)
            [x_sort,sort_ind] = sort(x(iexp,:));
            plot(x_sort,d(iexp,sort_ind),'k-')
        end
%         plot(d','k-')
        errorbar(mean(sort(x,2)),mean(d,1),ste(d,1),'.-','MarkerSize',20,'LineWidth',1)
        title(sprintf('%s Behavior-Response',av_label{iav+1}))
        figXAxis([],'n PCs',[0 max(x(:))+1],0:5:max(x(:)),0:5:max(x(:)))
        figYAxis([],'MI (bits)',[mi_lim_R])
        figAxForm
        subplot(2,3,3+(iav*3)); hold on
        d  = y_SB(:,:,iav+1);
        for iexp = 1:length(data)
            [x_sort,sort_ind] = sort(x(iexp,:));
            plot(x_sort,d(iexp,sort_ind),'k-')
        end
%         plot(d','k-')
        errorbar(mean(sort(x,2)),mean(d,1),ste(d,1),'.-','MarkerSize',20,'LineWidth',1)
        title(sprintf('%s Stimulus-Behavior',av_label{iav+1}))
        figXAxis([],'n PCs',[0 max(x(:))+1],0:5:max(x(:)),0:5:max(x(:)))
        figYAxis([],'MI (bits)',[mi_lim_SB])
        figAxForm        
    end
    if doLoadPrevious
        saveStr = ['MI_results_' runTag '_' prevTimestampID];
    else
        saveStr = ['MI_results_' runTag '_' timestampID];
    end
    print(fullfile(fnout,saveStr),'-dpdf','-fillpage')
end
