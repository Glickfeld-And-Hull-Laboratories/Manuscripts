mi_lim = [0 1.5];
ent_S_lim = [0.5 2.5];
ent_R_lim = [4 5.5];
attendFaceColor = [1,1,1;0,0,0];
fig_pid = figure;
suptitle('black=attn,white=no attn,color=mouse')
fig_mi = figure;
suptitle('black=attn,white=no attn,color=mouse')
fig_ent = figure;
suptitle('black=attn,white=no attn,color=mouse')

ms_expt = cellfun(@(x) x(1:3),{data.expt_name},'unif',0);
for iattend = 0:1
    if iattend == 0
        h = [];
        ms_label = [];
    end
    idx = cell2mat({data.hasAttention}) == iattend;
    figure(fig_mi)
    subplot 131; hold on
    mice = unique(ms_expt(idx)); 
    ms_label = cat(2,ms_label,mice);
    for im = 1:length(mice)
        idx_ms = strcmp(ms_expt,mice(im)) & idx;
        x = mean(PID.MI_SR.Vis(:,idx_ms),1);
        [xerr_l,xerr_u] = ciFromBoot(PID.MI_SR.Vis(:,idx_ms),95);
        y = mean(PID.MI_SR.Aud(:,idx_ms),1);
        [yerr_l,yerr_u] = ciFromBoot(PID.MI_SR.Aud(:,idx_ms),95);
        eb = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'o','MarkerSize',10);
        eb.MarkerFaceColor = attendFaceColor(iattend+1,:);
    end
    %     plot(x,y,'.','MarkerSize',ssize)
    title('mutual info in the stimulus and response')
    subplot 132; hold on
    for im = 1:length(mice)
        idx_ms = strcmp(ms_expt,mice(im)) & idx;
        x = mean(PID.MI_BR.Vis(:,idx_ms),1);
        [xerr_l,xerr_u] = ciFromBoot(PID.MI_BR.Vis(:,idx_ms),95);
        y = mean(PID.MI_BR.Aud(:,idx_ms),1);
        [yerr_l,yerr_u] = ciFromBoot(PID.MI_BR.Aud(:,idx_ms),95);
        eb = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'o','MarkerSize',10);
        eb.MarkerFaceColor = attendFaceColor(iattend+1,:);
    end
%     h(iattend+1) = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'.','MarkerSize',ssize);
    if iattend == 1
%         legend(attn_label,'location','northeast')
%         refline(1)
        for iplot = 1:3
        subplot(1,3,iplot)
        figXAxis([],'Visual',mi_lim)
        figYAxis([],'Auditory',mi_lim)
        figAxForm
        refline(1)
        end
    end
    title('mutual info in the choice and response')
    subplot 133; hold on
    for im = 1:length(mice)
        idx_ms = strcmp(ms_expt,mice(im)) & idx;
        x = mean(PID.MI_SB.Vis(:,idx_ms),1);
        [xerr_l,xerr_u] = ciFromBoot(PID.MI_SB.Vis(:,idx_ms),95);
        y = mean(PID.MI_SB.Aud(:,idx_ms),1);
        [yerr_l,yerr_u] = ciFromBoot(PID.MI_SB.Aud(:,idx_ms),95);
        eb = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'o','MarkerSize',10);
        eb.MarkerFaceColor = attendFaceColor(iattend+1,:);
    end
    title('mutual info in the stimulus and choice')
    
    
    figure(fig_pid)
    subplot 321; hold on
    for im = 1:length(mice)
        idx_ms = strcmp(ms_expt,mice(im)) & idx;
        x = mean(PID.II.Vis(:,idx_ms)./PID.MI_SR.Vis(:,idx_ms),1);
        [xerr_l,xerr_u] = ciFromBoot(PID.II.Vis(:,idx_ms)./PID.MI_SR.Vis(:,idx_ms),95);
        y = mean(PID.II.Aud(:,idx_ms)./PID.MI_SR.Aud(:,idx_ms),1);
        [yerr_l,yerr_u] = ciFromBoot(PID.II.Aud(:,idx_ms)./PID.MI_SR.Aud(:,idx_ms),95);
        eb = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'o','MarkerSize',10);
        eb.MarkerFaceColor = attendFaceColor(iattend+1,:);
    end
    axis([0 1 0 1])
    title('Percent Sensory info in Pop that drives behavior  (PID/SR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm
    if iattend == 1
        refline(1)
    end

    subplot 322; hold on
    for im = 1:length(mice)
        idx_ms = strcmp(ms_expt,mice(im)) & idx;
        x = mean(1-PID.II.Vis(:,idx_ms)./PID.MI_SR.Vis(:,idx_ms),1);
        [xerr_l,xerr_u] = ciFromBoot(1-PID.II.Vis(:,idx_ms)./PID.MI_SR.Vis(:,idx_ms),95);
        y = mean(1-PID.II.Aud(:,idx_ms)./PID.MI_SR.Aud(:,idx_ms),1);
        [yerr_l,yerr_u] = ciFromBoot(1-PID.II.Aud(:,idx_ms)./PID.MI_SR.Aud(:,idx_ms),95);
        eb = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'o','MarkerSize',10);
        eb.MarkerFaceColor = attendFaceColor(iattend+1,:);
    end
     axis([0 1 0 1])
    refline(1)
    title('Percent Sensory info in Pop that does not drives behavior (1-PID/SR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm

    subplot 323; hold on
    for im = 1:length(mice)
        idx_ms = strcmp(ms_expt,mice(im)) & idx;
        x = mean(1-PID.II.Vis(:,idx_ms)./PID.MI_BR.Vis(:,idx_ms),1);
        [xerr_l,xerr_u] = ciFromBoot(1-PID.II.Vis(:,idx_ms)./PID.MI_BR.Vis(:,idx_ms),95);
        y = mean(1-PID.II.Aud(:,idx_ms)./PID.MI_BR.Aud(:,idx_ms),1);
        [yerr_l,yerr_u] = ciFromBoot(1-PID.II.Aud(:,idx_ms)./PID.MI_BR.Aud(:,idx_ms),95);
        eb = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'o','MarkerSize',10);
        eb.MarkerFaceColor = attendFaceColor(iattend+1,:);
    end
    axis([0 1 0 1])
    refline(1)
    title('Percent Choice info in Pop not related to Stim (1-PID/BR)')
    xlabel('Visual Condition')
    ylabel('Auditory Condition')
    figAxForm

    subplot 324; hold on
    for im = 1:length(mice)
        idx_ms = strcmp(ms_expt,mice(im)) & idx;
        x = mean(1-PID.II.Vis(:,idx_ms)./PID.MI_SB.Vis(:,idx_ms),1);
        [xerr_l,xerr_u] = ciFromBoot(1-PID.II.Vis(:,idx_ms)./PID.MI_SB.Vis(:,idx_ms),95);
        y = mean(1-PID.II.Aud(:,idx_ms)./PID.MI_SB.Aud(:,idx_ms),1);
        [yerr_l,yerr_u] = ciFromBoot(1-PID.II.Aud(:,idx_ms)./PID.MI_SB.Aud(:,idx_ms),95);
        eb = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'o','MarkerSize',10);
        eb.MarkerFaceColor = attendFaceColor(iattend+1,:);
    end
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
    
    figure(fig_ent)
    for im = 1:length(mice)
        idx_ms = strcmp(ms_expt,mice(im)) & idx;
        x = mean(ENT.Vis_S(:,idx_ms),1);
        [xerr_l,xerr_u] = ciFromBoot(ENT.Vis_S(:,idx_ms),95);
        subplot 221; hold on
        y = mean(ENT.Vis_R(:,idx_ms),1);
        [yerr_l,yerr_u] = ciFromBoot(ENT.Vis_R(:,idx_ms),95);
        eb = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'o','MarkerSize',10);
        eb.MarkerFaceColor = attendFaceColor(iattend+1,:);
        subplot 223; hold on
        y = mean(PID.II.Vis(:,idx_ms)./PID.MI_SR.Vis(:,idx_ms),1);
        [yerr_l,yerr_u] = ciFromBoot(PID.II.Vis(:,idx_ms)./PID.MI_SR.Vis(:,idx_ms),95);
        eb = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'o','MarkerSize',10);
        eb.MarkerFaceColor = attendFaceColor(iattend+1,:);
        
        x = mean(ENT.Aud_S(:,idx_ms),1);
        [xerr_l,xerr_u] = ciFromBoot(ENT.Aud_S(:,idx_ms),95);
        subplot 222; hold on
        y = mean(ENT.Aud_R(:,idx_ms),1);
        [yerr_l,yerr_u] = ciFromBoot(ENT.Aud_R(:,idx_ms),95);
        eb = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'o','MarkerSize',10);
        eb.MarkerFaceColor = attendFaceColor(iattend+1,:);
        subplot 224; hold on
        y = mean(PID.II.Aud(:,idx_ms)./PID.MI_SR.Aud(:,idx_ms),1);
        [yerr_l,yerr_u] = ciFromBoot(PID.II.Aud(:,idx_ms)./PID.MI_SR.Aud(:,idx_ms),95);
        eb = errorbar(x,y,yerr_l,yerr_u,xerr_l,xerr_u,'o','MarkerSize',10);
        eb.MarkerFaceColor = attendFaceColor(iattend+1,:);
    end
    
    if iattend == 1
        x = mean(ENT.Vis_S,1);
        y = mean(ENT.Vis_R,1);
        [r,p] = corrcoef(x,y);
        r_ent = r(1,2);
        p_ent = p(1,2);
        x = mean(ENT.Aud_S,1);
        y = mean(ENT.Aud_R,1);
        [r,p] = corrcoef(x,y);
        r_ent = cat(2,r_ent,r(1,2));
        p_ent = cat(2,p_ent,p(1,2));
        
        x = mean(ENT.Vis_S,1);
        y = mean(PID.II.Vis./PID.MI_SR.Vis,1);
        [r,p] = corrcoef(x,y);
        r_pid = r(1,2);
        p_pid = p(1,2);
        x = mean(ENT.Aud_S,1);
        y = mean(PID.II.Aud./PID.MI_SR.Aud,1);
        [r,p] = corrcoef(x,y);
        r_pid = cat(2,r_pid,r(1,2));
        p_pid = cat(2,p_pid,p(1,2));
        
        for iav = 1:2
            subplot(2,2,iav)
            figXAxis([],'Entropy in Stim',ent_S_lim)
            figYAxis([],'Entropy in Resp',ent_R_lim)
            figAxForm
            title(sprintf('%s trials, r=%s,p=%s',...
                av_label{iav},sigfigString(r_ent(iav)),...
                sigfigString(p_ent(iav))))   
            
            
            subplot(2,2,iav+2)
            figXAxis([],'Entropy in Stim',ent_S_lim)
            figYAxis([],'PID/MI_SR',[0 1])
            figAxForm
            title(sprintf('%s trials, r=%s,p=%s',...
                av_label{iav},sigfigString(r_pid(iav)),...
                sigfigString(p_pid(iav))))
        end
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
figure(fig_ent)
print(fullfile(fnout,['ENT' saveStr]),'-dpdf','-fillpage')
