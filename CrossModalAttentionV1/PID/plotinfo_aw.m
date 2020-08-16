% attend = 1;
cc(1,:)=[0 0 0];
cc(2,:)=[0.75 0 0.75];
attn_label = {'no attn';'attn'};
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
    ENT = struct;
    ENT.Vis_S = nan(nBoot,length(D));
    ENT.Vis_B = nan(nBoot,length(D));
    ENT.Vis_R = nan(nBoot,length(D));
    ENT.Aud_S = nan(nBoot,length(D));
    ENT.Aud_B = nan(nBoot,length(D));
    ENT.Aud_R = nan(nBoot,length(D));
    for iboot = 1:nBoot
        PID.II.Vis(iboot,:) = PIDresults.PID{iboot}.II.Vis;
        PID.MI_SR.Vis(iboot,:) = PIDresults.PID{iboot}.MI_SR.Vis;
        PID.MI_BR.Vis(iboot,:) = PIDresults.PID{iboot}.MI_BR.Vis;
        PID.MI_SB.Vis(iboot,:) = PIDresults.PID{iboot}.MI_SB.Vis;
        PID.II.Aud(iboot,:) = PIDresults.PID{iboot}.II.Aud;
        PID.MI_SR.Aud(iboot,:) = PIDresults.PID{iboot}.MI_SR.Aud;
        PID.MI_BR.Aud(iboot,:) = PIDresults.PID{iboot}.MI_BR.Aud;
        PID.MI_SB.Aud(iboot,:) = PIDresults.PID{iboot}.MI_SB.Aud;  

        ENT.Vis_S(iboot,:) = PIDresults.ENT{iboot}.Vis_S;
        ENT.Vis_B(iboot,:) = PIDresults.ENT{iboot}.Vis_B;
        ENT.Vis_R(iboot,:) = PIDresults.ENT{iboot}.Vis_R;
        ENT.Aud_S(iboot,:) = PIDresults.ENT{iboot}.Aud_S;
        ENT.Aud_B(iboot,:) = PIDresults.ENT{iboot}.Aud_B;
        ENT.Aud_R(iboot,:) = PIDresults.ENT{iboot}.Aud_R;
    end
else
    PID = PIDresults.PID;
end
%% plot PID results across experiments
if ~doMultiPC & ~doMultiCluster
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
        h_scat = [];
    end
    idx = cell2mat({data.hasAttention}) == iattend;
    figure(fig_mi)
    subplot 131; hold on
    plot(PID.MI_SR.Vis(idx),PID.MI_SR.Aud(idx),'.','MarkerSize',ssize)
    title('mutual info in the stimulus and response')
    subplot 132; hold on
    plot(PID.MI_BR.Vis(idx),PID.MI_BR.Aud(idx),'.','MarkerSize',ssize)
    title('mutual info in the choice and response')
    subplot 133; hold on
    h(iattend+1) = plot(PID.MI_BR.Vis(idx),PID.MI_BR.Aud(idx),'.','MarkerSize',ssize);
    title('mutual info in the stimulus and choice')
    if iattend == 1
        legend(h,attn_label)
    end
    for iplot = 1:3
        subplot(1,3,iplot)
        figXAxis([],'Visual',mi_lim)
        figYAxis([],'Auditory',mi_lim)
        figAxForm
        refline(1)
    end
    
    figure(fig_pid)
    subplot 321; hold on
    h_scat(iattend+1) = scatter(...
        PID.II.Vis(idx)./PID.MI_SR.Vis(idx),PID.II.Aud(idx)./PID.MI_SR.Aud(idx),ssize);
    axis([0 1 0 1])
    title('Percent Sensory info in Pop that drives behavior  (PID/SR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm
    if iattend == 1
        legend(h_scat,attn_label,'location','southwestoutside')
        refline(1)
    end
    if iattend == 0
        idx_halfattend = cellfun(@(x) strcmp(x(1:3),'750'),{data.expt_name});
        plot(PID.II.Vis(idx_halfattend)./PID.MI_SR.Vis(idx_halfattend),...
            PID.II.Aud(idx_halfattend)./PID.MI_SR.Aud(idx_halfattend),'k.','MarkerSize',5)
    end

    subplot 322; hold on
    scatter(1-PID.II.Vis(idx)./PID.MI_SR.Vis(idx),1-PID.II.Aud(idx)./PID.MI_SR.Aud(idx),ssize)
    axis([0 1 0 1])
    refline(1)
    title('Percent Sensory info in Pop that does not drives behavior (1-PID/SR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm
    if iattend == 0
        plot(1-PID.II.Vis(idx_halfattend)./PID.MI_SR.Vis(idx_halfattend),1-PID.II.Aud(idx_halfattend)./PID.MI_SR.Aud(idx_halfattend),'k.','MarkerSize',5)
    end

    subplot 323; hold on
    scatter(1-PID.II.Vis(idx)./PID.MI_BR.Vis(idx),1-PID.II.Aud(idx)./PID.MI_BR.Aud(idx),ssize)
    axis([0 1 0 1])
    refline(1)
    title('Percent Choice info in Pop not related to Stim (1-PID/BR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm
    if iattend == 0
        plot(1-PID.II.Vis(idx_halfattend)./PID.MI_BR.Vis(idx_halfattend),1-PID.II.Aud(idx_halfattend)./PID.MI_BR.Aud(idx_halfattend),'k.','MarkerSize',5)
    end

    subplot 324; hold on
    scatter(1-PID.II.Vis(idx)./PID.MI_SB.Vis(idx),1-PID.II.Aud(idx)./PID.MI_SB.Aud(idx),ssize)
    axis([0 1 0 1])
    refline(1)
    title('Percent Stimulus info that drives behavior via an alt. pathway (1-PID/SB)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm
    if iattend == 0
        plot(1-PID.II.Vis(idx_halfattend)./PID.MI_SB.Vis(idx_halfattend),1-PID.II.Aud(idx_halfattend)./PID.MI_SB.Aud(idx_halfattend),'k.','MarkerSize',5)
    end

%     figure(fig_m_pid)
    subplot(3,2,iattend+5); hold on

    % non_readout_sensory_info=S_R_info-I_II;
    b1=mean([(PID.II.Vis(idx)./PID.MI_SR.Vis(idx))',(PID.II.Aud(idx)./PID.MI_SR.Aud(idx))']);
    b1err = ste([(PID.II.Vis(idx)./PID.MI_SR.Vis(idx))',(PID.II.Aud(idx)./PID.MI_SR.Aud(idx))'],1);
    b1=1-b1;% sensory information in neural response R that is not read out for behavior

    %internal_choice_info=C_R_info-I_II;
    b2=mean([(PID.II.Vis(idx)./PID.MI_BR.Vis(idx))',(PID.II.Aud(idx)./PID.MI_BR.Aud(idx))']);
    b2err = ste([(PID.II.Vis(idx)./PID.MI_BR.Vis(idx))',(PID.II.Aud(idx)./PID.MI_BR.Aud(idx))'],1);
    b2=1-b2;% choice information in neural response R that is not related to the stimulus

    %S_C_info_from_unobserved_R=S_C_info-I_II;
    b3=mean([(PID.II.Vis(idx)./PID.MI_SB.Vis(idx))',(PID.II.Aud(idx)./PID.MI_SB.Aud(idx))']);
    b3err = ste([(PID.II.Vis(idx)./PID.MI_SB.Vis(idx))',(PID.II.Aud(idx)./PID.MI_SB.Aud(idx))'],1);
    b3 = 1-b3; % the part of "behavioral performance" that cannot be explained 
               % with recorded neural feature R
        h_sum = [];
    for iav = 1:2
        h_sum{iav} = errorbar([b1(iav);b2(iav);b3(iav);],[b1err(iav);b2err(iav);b3err(iav);],'.','MarkerSize',20);
        h_sum{iav}.Color = cc(iav,:);
    end
    if iattend == 1
        legend([h_sum{1},h_sum{2}],av_label,'location','northeastoutside')
    end
    figXAxis([],'',[0 4],1:3,{'1-II/SR';'1-II/BR';'1-II/SB'})
    xtickangle(-45)
    figYAxis([],'II',[0 1])
    figAxForm
    title({attn_label{iattend+1};'error across expt'})
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

%% compare npcs or nclusters used
doExptPlot = false;
% exptNPCsID = '200721_1946';
mi_lim_R = [0 1.3];
mi_lim_SB = [0 0.75];
if doMultiPC | doMultiCluster
    if doMultiPC
        n = npcs;
    elseif doMultiCluster
        n = nc;
    end
%     load(fullfile(fnout,['PID_results_' exptNPCsID]));
    y_SR = nan(length(data),length(n),2);
    y_BR = nan(length(data),length(n),2);
    y_SB = nan(length(data),length(n),2);
    x = n;
    for iexp = 1:length(data)
        if doExptPlot
            figure
            suptitle(data(iexp).expt_name)
        end
        for iav = 0:1
            if iav == 0
                ntrials = length(PIDresults.D{1}{iexp}.Vis.Rlabels);
            else
                ntrials = length(PIDresults.D{1}{iexp}.Aud.Rlabels);
            end
            for i = 1:length(n)
                if iav == 0
                    y_SR(iexp,i,iav+1) = PIDresults.PID{i}.MI_SR.Vis(iexp);
                    y_BR(iexp,i,iav+1) = PIDresults.PID{i}.MI_BR.Vis(iexp);
                    y_SB(iexp,i,iav+1) = PIDresults.PID{i}.MI_SB.Vis(iexp);
                else
                    y_SR(iexp,i,iav+1) = PIDresults.PID{i}.MI_SR.Aud(iexp);
                    y_BR(iexp,i,iav+1) = PIDresults.PID{i}.MI_BR.Aud(iexp);
                    y_SB(iexp,i,iav+1) = PIDresults.PID{i}.MI_SB.Aud(iexp);
                end
            end
            if doExptPlot
                subplot(2,3,1+(iav*3))
                plot(y_SR(iexp,:,iav+1),'.-','MarkerSize',20,'LineWidth',1)
                title(sprintf('%s Stimulus-Response, %s trials',av_label{iav+1},num2str(ntrials)))
                figXAxis([],'n PCs',[0 length(n)+2],1:length(n)+1,[n,0])
                figYAxis([],'MI (bits)',[])
                figAxForm
                subplot(2,3,2+(iav*3))
                plot(y_BR(iexp,:,iav+1),'.-','MarkerSize',20,'LineWidth',1)
                title(sprintf('%s Behavior-Response, %s trials',av_label{iav+1},num2str(ntrials)))
                figXAxis([],'n PCs',[0 length(n)+2],1:length(n)+1,[n,0])
                figYAxis([],'MI (bits)',[])
                figAxForm
                subplot(2,3,3+(iav*3))
                plot(y_SB(iexp,:,iav+1),'.-','MarkerSize',20,'LineWidth',1)
                title(sprintf('%s Stimulus-Behavior, %s trials',av_label{iav+1},num2str(ntrials)))
                figXAxis([],'n PCs',[0 length(n)+2],1:length(n)+1,[n,0])
                figYAxis([],'MI (bits)',[])
                figAxForm
            end
        end
    end
    figure
    for iav = 0:1
        subplot(2,3,1+(iav*3)); hold on
        d  = y_SR(:,:,iav+1);
        for iexp = 1:length(data)
%             [x_sort,sort_ind] = sort(x(iexp,:));
            plot(x,d(iexp,:),'k-')
        end
%         plot(d','k-')
        errorbar(x,mean(d,1),ste(d,1),'.-','MarkerSize',20,'LineWidth',1)
        if doMultiPC
            figXAxis([],'n PCs',[0 max(x(:))+1],0:5:max(x(:)),0:5:max(x(:)))
        elseif doMultiCluster
            figXAxis([],'n Clusters',[0 max(x(:))+1],x,x)
        end
        figYAxis([],'MI (bits)',[mi_lim_R])
        figAxForm
        title(sprintf('%s Stimulus-Response',av_label{iav+1}))
        
        subplot(2,3,2+(iav*3)); hold on
        d  = y_BR(:,:,iav+1);        
        for iexp = 1:length(data)
%             [x_sort,sort_ind] = sort(x(iexp,:));
            plot(x,d(iexp,:),'k-')
        end
%         plot(d','k-')
        errorbar(x,mean(d,1),ste(d,1),'.-','MarkerSize',20,'LineWidth',1)
        if doMultiPC
            figXAxis([],'n PCs',[0 max(x(:))+1],0:5:max(x(:)),0:5:max(x(:)))
        elseif doMultiCluster
            figXAxis([],'n Clusters',[0 max(x(:))+1],x,x)
        end
        figYAxis([],'MI (bits)',[mi_lim_R])
        figAxForm
        title(sprintf('%s Behavior-Response',av_label{iav+1}))
        
        subplot(2,3,3+(iav*3)); hold on
        d  = y_SB(:,:,iav+1);
        for iexp = 1:length(data)
%             [x_sort,sort_ind] = sort(x(iexp,:));
            plot(x,d(iexp,:),'k-')
        end
%         plot(d','k-')
        errorbar(x,mean(d,1),ste(d,1),'.-','MarkerSize',20,'LineWidth',1)
        if doMultiPC
            figXAxis([],'n PCs',[0 max(x(:))+1],0:5:max(x(:)),0:5:max(x(:)))
        elseif doMultiCluster
            figXAxis([],'n Clusters',[0 max(x(:))+1],x,x)
        end
        figYAxis([],'MI (bits)',[mi_lim_R])
        figAxForm   
        title(sprintf('%s Stimulus-Behavior',av_label{iav+1}))     
    end
    if doLoadPrevious
        saveStr = ['MI_results_' runTag '_' prevTimestampID];
    else
        saveStr = ['MI_results_' runTag '_' timestampID];
    end
    print(fullfile(fnout,saveStr),'-dpdf','-fillpage')
end
%% stimulus entropy
if ~doMultiPC | ~doMultiCluster | ~doBootstrapPID
if binaryS
    ent_S_lim = [0.5 1];
    ent_R_lim = [1 3.2];
else
    ent_S_lim = [0.5 2];
    ent_R_lim = [3 6];
end
figure
for iattend = 0:1
    idx = cell2mat({data.hasAttention}) == iattend;

    subplot 321; hold on
    x = PIDresults.ENT.Vis_S(idx);
    y = PIDresults.ENT.Vis_R(idx);
    plot(x,y,'.','MarkerSize',20);
    subplot 323; hold on
    y = PIDresults.PID.MI_SR.Vis(idx);
    plot(x,y,'.','MarkerSize',20);
    subplot 325; hold on
    y = PIDresults.PID.II.Vis(idx);
    plot(x,y,'.','MarkerSize',20);
    
    subplot 322; hold on
    x = PIDresults.ENT.Aud_S(idx);
    y = PIDresults.ENT.Aud_R(idx);
    plot(x,y,'.','MarkerSize',20);
    subplot 324; hold on
    y = PIDresults.PID.MI_SR.Aud(idx);
    plot(x,y,'.','MarkerSize',20);
    subplot 326; hold on
    y = PIDresults.PID.II.Aud(idx);
    plot(x,y,'.','MarkerSize',20);
    
end
for iav = 1:2
    subplot(3,2,iav)
    refline(1)
    figXAxis([],'Entropy in Stim',ent_S_lim)
    figYAxis([],'Entropy in Resp',ent_R_lim)
    figAxForm
    title(av_label{iav})
    
    subplot(3,2,iav+2)
    figXAxis([],'Entropy in Stim',ent_S_lim)
    figYAxis([],'MI (bits) SR',[])
    figAxForm
    title(av_label{iav})
    
    subplot(3,2,iav+4)
    figXAxis([],'Entropy in Stim',ent_S_lim)
    figYAxis([],'PID',[])
    figAxForm
    title(av_label{iav})
end
    if doLoadPrevious
        saveStr = ['StimEntropyCorrs_' runTag '_' prevTimestampID];
    else
        saveStr = ['StimEntropyCorr_' runTag '_' timestampID];
    end
    print(fullfile(fnout,saveStr),'-dpdf','-fillpage')
end