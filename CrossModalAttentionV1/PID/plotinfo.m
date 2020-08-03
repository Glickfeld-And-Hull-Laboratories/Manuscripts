
cc(1,:)=[0 0 1];  % no Attention is BLUE
cc(2,:)=[1 0 0];  % Attention is RED
number=false;
ssize=20*ones(size(attend));

figure
subplot(2,2,1), scatter(PID.II.Vis./PID.MI_SR.Vis,PID.II.Aud./PID.MI_SR.Aud,ssize,cc(attend+1,:))
axis([0 1 0 1])
refline(1)
title('Percent Sensory info in Pop that drives behavior (Attend==red)')
xlabel('Visual Condition')
ylabel('Auditory Condition')

subplot(2,2,2), scatter(1-PID.II.Vis./PID.MI_SR.Vis,1-PID.II.Aud./PID.MI_SR.Aud,ssize,cc(attend+1,:))
axis([0 1 0 1])
refline(1)
title('Percent Sensory info in Pop that does not drives behavior')
xlabel('Visual Condition')
ylabel('Auditory Condition')

subplot(2,2,3), scatter(1-PID.II.Vis./PID.MI_BR.Vis,1-PID.II.Aud./PID.MI_BR.Aud,ssize,cc(attend+1,:))
axis([0 1 0 1])
refline(1)
title('Percent Choice info in Pop not related to Stim')
xlabel('Visual Condition')
ylabel('Auditory Condition')

subplot(2,2,4), scatter(1-PID.II.Vis./PID.MI_SB.Vis,1-PID.II.Aud./PID.MI_SB.Aud,ssize,cc(attend+1,:))
axis([0 1 0 1])
refline(1)
title('Percent Stimulus info that drives behavior via an alt. pathway')
xlabel('Visual Condition')
ylabel('Auditory Condition')

figure

% non_readout_sensory_info=S_R_info-I_II;
b1=mean([(PID.II.Vis./PID.MI_SR.Vis)',(PID.II.Aud./PID.MI_SR.Aud)']);
b1=1-b1;% sensory information in neural response R that is not read out for behavior

%internal_choice_info=C_R_info-I_II;
b2=mean([(PID.II.Vis./PID.MI_BR.Vis)',(PID.II.Aud./PID.MI_BR.Aud)']);
b2=1-b2;% choice information in neural response R that is not related to the stimulus

%S_C_info_from_unobserved_R=S_C_info-I_II;
b3=mean([(PID.II.Vis./PID.MI_SB.Vis)',(PID.II.Aud./PID.MI_SB.Aud)']);
b3 = 1-b3; % the part of "behavioral performance" that cannot be explained 
           % with recorded neural feature R

bar([b1;b2;b3;]) 

figure
hist([ENT.Aud_S(attend==0)',ENT.Vis_S(attend==0)'])
xlabel('Stimulus Entropy')
ylabel('count')
legend('Auditory Trials','Visual Trials')
title('Un Attending')
figure
hist([ENT.Aud_S(attend==1)',ENT.Vis_S(attend==1)'])
title('Attending')
xlabel('Stimulus Entropy')
ylabel('count')
legend('Auditory Trials','Visual Trials')

% figure
% h=subplot(2,2,1), scatter(PID.MI_SB.Vis,PID.MI_SB.Aud,ssize,cc(attend+1,:))
% h1=get(h,'Children');
% x=get(h1,'Xdata');
% y=get(h1,'Ydata');
% a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
% if(number)
%     text(h,x, y, c);
%     h1.CData=[1 1 1];
% end
% m=max(1,max([x,y]));
% axis([0 m 0 m])
% refline(1)
% title('Stim info in behavior (Attend==red)')
% xlabel('Visual Information in Behavior')
% ylabel('Auditory Information in Behavior')
% 
% 
% h=subplot(2,2,2), scatter(V_R_info./V_B_info,A_R_info./A_B_info,ssize,cc(attend+1,:));
% h1=get(h,'Children');
% x=get(h1,'Xdata');
% y=get(h1,'Ydata');
% a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
% if(number)
%     text(h,x, y, c);
%     h1.CData=[1 1 1];
% end
% m=max(1,max([x,y]));
% axis([0 m 0 m])
% refline(1)
% xlabel('Percent Visual Cue Info')
% ylabel('Percent Auditory Cue Info')
% 
% h=subplot(2,2,3), scatter(Vinternal_choice_info./VB_R_info,Ainternal_choice_info./AB_R_info,ssize,cc(attend+1,:));
% h1=get(h,'Children');
% x=get(h1,'Xdata');
% y=get(h1,'Ydata');
% a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
% if(number)
%     text(h,x, y, c);
%     h1.CData=[1 1 1];
% end
% m=max(1,max([x,y]));
% axis([0 m 0 m])
% refline(1)
% xlabel('Percent Pop resp that drives chioce independent of Visual stim')
% ylabel('Percent Pop resp that drives chioce independent of Auditory stim')
% title('Visual information correlates more with behavior')
% 
% h=subplot(2,2,4), scatter(Vnon_readout_sensory_info./V_R_info,Anon_readout_sensory_info./A_R_info,ssize,cc(attend+1,:));
% h1=get(h,'Children');
% x=get(h1,'Xdata');
% y=get(h1,'Ydata');
% a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
% if(number)
%     text(h,x, y, c);
%     h1.CData=[1 1 1];
% end
% m=max(1,max([x,y]));
% axis([0 m 0 m])
% refline(1)
% xlabel('Non-Readout Visual Info over total Visual Information')
% ylabel('Non-Readout Auditory Info over total Auditory Information')
% title('Alot of the Auditory information is not used to drive behavior')
% 
% 
% figure
% h=subplot(2,2,1), scatter(V_B_info,A_B_info,ssize,cc(attend+1,:))
% h1=get(h,'Children');
% x=get(h1,'Xdata');
% y=get(h1,'Ydata');
% a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
% if(number)
%     text(h,x, y, c);
%     h1.CData=[1 1 1];
% end
% m=max(max([x,y]));
% axis([0 m 0 m])
% refline(1)
% title('Stim info in behavior (0.4175 bits = 75 % correct)')
% xlabel('Visual Information in Behavior')
% ylabel('Auditory Information in Behavior')
% legend('Unattended','Attended')
% 
% 
% h=subplot(2,2,2), scatter(V_R_info,A_R_info,ssize,cc(attend+1,:));
% h1=get(h,'Children');
% x=get(h1,'Xdata');
% y=get(h1,'Ydata');
% a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
% if(number)
%     text(h,x, y, c);
%     h1.CData=[1 1 1];
% end
% m=max(max([x,y]));
% axis([0 m 0 m])
% refline(1)
% xlabel('Bits Visual Cue Info')
% ylabel('Bits Auditory Cue Info')
% title('There is more visual info than auditory info in the population response')
% 
% h=subplot(2,2,3), scatter(Vinternal_choice_info,Ainternal_choice_info,ssize,cc(attend+1,:));
% h1=get(h,'Children');
% x=get(h1,'Xdata');
% y=get(h1,'Ydata');
% a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
% if(number)
%     text(h,x, y, c);
%     h1.CData=[1 1 1];
% end
% m=max(max([x,y]));
% axis([0 m 0 m])
% refline(1)
% xlabel('Bits Pop resp that drives chioce independent of Visual stim')
% ylabel('Bits Pop resp that drives chioce independent of Auditory stim')
% title('Visual information correlates more with behavior')
% 
% h=subplot(2,2,4), scatter(Vnon_readout_sensory_info,Anon_readout_sensory_info,ssize,cc(attend+1,:));
% h1=get(h,'Children');
% x=get(h1,'Xdata');
% y=get(h1,'Ydata');
% a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
% if(number)
%     text(h,x, y, c);
%     h1.CData=[1 1 1];
% end
% m=max(max([x,y]));
% axis([0 m 0 m])
% refline(1)
% xlabel('Non-Readout Visual Info')
% ylabel('Non-Readout Auditory Info')
% title('Alot of the Auditory information is not used to drive behavior')
% 
% 
% 
