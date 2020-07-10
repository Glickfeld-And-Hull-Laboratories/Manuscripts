%%% ATTENDED VS UNATTENDED 
%%% trial type cued by presence or absense of tone.  
%%% invalid visual trial means auditory target comes
%%% invalid auditory trial means visual comes first
%%% 
%%% a mouse has attention if an invalid stimulus does not elilct a response
%%% (as often).    These were rare.  

clear all

load('attnData_neuronResp_trialOutcome.mat')
addpath ../PID
addpath ../PID/glpkmex/
addpath ../PID/glpkmex/win64/

%Add frame number as a categorical regressor for predicting behavior.

for n=1:length(data);
attend(n)=data(n).hasAttention;
NT=length(data(n).dff);

D=size(data(n).dff{1},2);

alldata=[];
resp=[];
vstim=[];
astim=[];
frame=[];
cue=[];
for t=1:NT
    alldata=[alldata;data(n).dff{t}];
    np=size(data(n).dff{t},1);
    if(data(n).visualTrial(t)) % then visual trial
        if(data(n).mouseResponse{t}=='m')
            resp=[resp;zeros(np,1);];
            vstim=[vstim;zeros(np-1,1);data(n).stimulus_orientation(t)];
        elseif(data(n).mouseResponse{t}=='fa')
            resp=[resp;zeros(np-1,1);1;];        
            vstim=[vstim;zeros(np,1)];
        elseif(data(n).mouseResponse{t}=='h')
            resp=[resp;zeros(np-1,1);1;];
            vstim=[vstim;zeros(np-1,1);data(n).stimulus_orientation(t)];
        end
        astim=[astim;zeros(np,1)];
        frame=[frame;[1:np]';];
        cue=[cue;ones(np,1);];% for visual trials
        
    else %auditory trial
        if(data(n).mouseResponse{t}=='m')
            resp=[resp;zeros(np,1);];
            astim=[astim;zeros(np-1,1);1];
        elseif(data(n).mouseResponse{t}=='fa')
            resp=[resp;zeros(np-1,1);1;];        
            astim=[astim;zeros(np,1)];
        elseif(data(n).mouseResponse{t}=='h')
            resp=[resp;zeros(np-1,1);1;];
            astim=[astim;zeros(np-1,1);data(n).stimulus_soundVolume(t)];
        end
        vstim=[vstim;zeros(np,1)];
        frame=[frame;[1:np]';]; 
        cue=[cue;zeros(np,1);];%0 for auditory trials

    end
        
end
[C,IA,IC]=unique(vstim);
vstim=IC-1;  % task difficulty 0 = distractor and 1 to 6 is 16 to 90 degrees.


astim(isnan(astim))=0;
[C,IA,IC]=unique(astim);
astim=IC-1; % no sound == 0;

% vstim(vstim>0)=1;
% astim(astim>0)=1;
vstim(vstim>5)=5;
idx = frame>3 & cue==1;% & vstim<=5; 
cue1=cue(idx);
vstim1=vstim(idx); 
astim1=astim(idx);
frame1=frame(idx);
resp1=resp(idx);
alldata1=alldata(idx,:);

astim(astim>5)=5;
idx = frame>3 & cue==0 ;%& astim<=5; %& cue==1;% & vstim<=3;
cue2=cue(idx);
vstim2=vstim(idx); 
astim2=astim(idx);
frame2=frame(idx)-1;
resp2=resp(idx);
alldata2=alldata(idx,:);

cue=[cue1;cue2];
vstim=[vstim1;vstim2];
astim=[astim1;astim2];
frame=[frame1;frame2];
resp=[resp1;resp2];
alldata=[alldata1;alldata2;];
npcs=15;
% 
%%%PCA reduction
[coeff,score]=pca(alldata);
score=score(:,1:min(npcs,size(alldata,2)));
alldata=zscore(score);
idx1=[ones(size(alldata1,1),1);zeros(size(alldata2,1),1);]==1;
idx2=~idx1;
alldata1=alldata(idx1,:);
alldata2=alldata(idx2,:);

% %Cluster reduction
% [coeff,score]=pca(alldata);
% score=score(:,1:min(npcs,size(alldata,2)));
% nc=50;
% [z,c]=kmeans(score,nc);
% datatemp=zeros(size(alldata,1),nc);
% for i=1:size(alldata,1)
%     datatemp(i,z(i))=1;
% end
% idx=sum(datatemp,1)>0;
% datatemp=datatemp(:,idx);
% idx1=[ones(size(alldata1,1),1);zeros(size(alldata2,1),1);]==1;
% idx2=~idx1;
% alldata=datatemp;
% alldata1=alldata(idx1,:);
% alldata2=alldata(idx2,:);
% idx=sum(alldata1,1)>1;
% alldata1=alldata1(:,idx);
% idx=sum(alldata2,1)>1;
% alldata2=alldata2(:,idx);

ns=length(unique(vstim1));
%[MI_r_vb,~,~,pvbgr] = calcMI(alldata1,(vstim1>0)+2*resp1+1);
[MI_r_vb,~,~,pvbgr] = calcMI(alldata1,vstim1+ns*resp1+1);
m=max(vstim1+ns*resp1+1);
if(m<2*ns)
    pvbgr(:,m+1:2*ns)=0;
end

ns=length(unique(astim2));
%[MI_r_ab,~,~,pabgr] = calcMI(alldata2,(astim2>0)+2*resp2+1);
[MI_r_ab,~,~,pabgr] = calcMI(alldata2,astim2+ns*resp2+1);
m=max(astim2+ns*resp2+1);
if(m<2*ns)
    pabgr(:,m+1:2*ns)=0;
end
    
% clear stim
% stim(vstim==0 & astim==0,1)=0;
% stim(vstim>0 & astim==0,1)=1;
% stim(vstim==0 & astim==1,1)=2;
% %stim(vstim>0 & astim==1,1)=3;
% 
% [MI_r_sb,~,~,p] = calcMI(alldata,stim+3*resp+1);
% 
% 
% %[MI_r_sb,~,~,p] = calcMI(alldata,stim+4*resp+1);


MI_rv(n) = calcMI(alldata1,vstim1+1);
MI_ra(n) = calcMI(alldata2,astim2+1);
MI_vb(n) = calcMI(resp1,vstim1+1);
MI_ab(n) = calcMI(resp2,astim2+1);
MI_Vrb(n) = calcMI(alldata1,resp1+1);
MI_Arb(n) = calcMI(alldata2,resp2+1);

[C,IA,IC]=unique(vstim1);
pIA=IA/sum(IA);
HV(n)=-sum(pIA.*log2(pIA));

[C,IA,IC]=unique(astim2);
pIA=IA/sum(IA);
HA(n)=-sum(pIA.*log2(pIA));

[C,IA,IC]=unique(resp1);
pIA=IA/sum(IA);
HBV(n)=-sum(pIA.*log2(pIA));

[C,IA,IC]=unique(resp2);
pIA=IA/sum(IA);
HBA(n)=-sum(pIA.*log2(pIA));



nc=50;
[z,c]=kmeans(pvbgr,nc);
%[z,c]=kmeans(alldata1,nc);
zc=hist(z,nc);
zc=zc/sum(zc);
ns=length(unique(vstim1));
pv=zeros(ns,nc,2);
for s=1:ns
for b=1:2
for r=1:nc    
%    pv(s,r,b)=zc(r)*c(r,s+ns*(b-1));    
    pv(s,r,b)=sum(pvbgr(z==r,s+ns*(b-1)));
end
end
end

[z,c]=kmeans(pabgr,nc);
%[z,c]=kmeans(alldata2,nc);
zc=hist(z,nc);
zc=zc/sum(zc);
ns=length(unique(astim2));
pa=zeros(ns,nc,2);
for s=1:ns
for b=1:2
for r=1:nc
    pa(s,r,b)=sum(pabgr(z==r,s+ns*(b-1)));
%    pa(s,r,b)=zc(r)*c(r,s+ns*(b-1));    
end
end
end

[VI_II(n), V_R_info(n), VB_R_info(n), V_B_info(n), Vnon_readout_sensory_info(n), Vinternal_choice_info(n), ...
V_B_info_from_unobserved_R(n)]=intersection_information(pv/sum(pv(:)));

[AI_II(n), A_R_info(n), AB_R_info(n), A_B_info(n), Anon_readout_sensory_info(n), Ainternal_choice_info(n), ...
A_B_info_from_unobserved_R(n)]=intersection_information(pa/sum(pa(:)));

[n,VI_II(n),AI_II(n)]

end
cc(1,:)=[0 0 1];
cc(2,:)=[1 0 0];
number=false;
ssize=20*ones(size(attend));
figure
h=subplot(2,2,1), scatter(V_B_info,A_B_info,ssize,cc(attend+1,:))
h1=get(h,'Children');
x=get(h1,'Xdata');
y=get(h1,'Ydata');
a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
if(number)
    text(h,x, y, c);
    h1.CData=[1 1 1];
end
m=max(1,max([x,y]));
axis([0 m 0 m])
refline(1)
title('Stim info in behavior (0.4175 bits = 75 % correct)')
xlabel('Visual Information in Behavior')
ylabel('Auditory Information in Behavior')
legend('Unattended','Attended')


h=subplot(2,2,2), scatter(V_R_info./V_B_info,A_R_info./A_B_info,ssize,cc(attend+1,:));
h1=get(h,'Children');
x=get(h1,'Xdata');
y=get(h1,'Ydata');
a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
if(number)
    text(h,x, y, c);
    h1.CData=[1 1 1];
end
m=max(1,max([x,y]));
axis([0 m 0 m])
refline(1)
xlabel('Percent Visual Cue Info')
ylabel('Percent Auditory Cue Info')
title('There is more visual info than auditory info in the population response')

h=subplot(2,2,3), scatter(Vinternal_choice_info./VB_R_info,Ainternal_choice_info./AB_R_info,ssize,cc(attend+1,:));
h1=get(h,'Children');
x=get(h1,'Xdata');
y=get(h1,'Ydata');
a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
if(number)
    text(h,x, y, c);
    h1.CData=[1 1 1];
end
m=max(1,max([x,y]));
axis([0 m 0 m])
refline(1)
xlabel('Percent Pop resp that drives chioce independent of Visual stim')
ylabel('Percent Pop resp that drives chioce independent of Auditory stim')
title('Visual information correlates more with behavior')

h=subplot(2,2,4), scatter(Vnon_readout_sensory_info./V_R_info,Anon_readout_sensory_info./A_R_info,ssize,cc(attend+1,:));
h1=get(h,'Children');
x=get(h1,'Xdata');
y=get(h1,'Ydata');
a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
if(number)
    text(h,x, y, c);
    h1.CData=[1 1 1];
end
m=max(1,max([x,y]));
axis([0 m 0 m])
refline(1)
xlabel('Non-Readout Visual Info over total Visual Information')
ylabel('Non-Readout Auditory Info over total Auditory Information')
title('Alot of the Auditory information is not used to drive behavior')


figure
h=subplot(2,2,1), scatter(V_B_info,A_B_info,ssize,cc(attend+1,:))
h1=get(h,'Children');
x=get(h1,'Xdata');
y=get(h1,'Ydata');
a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
if(number)
    text(h,x, y, c);
    h1.CData=[1 1 1];
end
m=max(max([x,y]));
axis([0 m 0 m])
refline(1)
title('Stim info in behavior (0.4175 bits = 75 % correct)')
xlabel('Visual Information in Behavior')
ylabel('Auditory Information in Behavior')
legend('Unattended','Attended')


h=subplot(2,2,2), scatter(V_R_info,A_R_info,ssize,cc(attend+1,:));
h1=get(h,'Children');
x=get(h1,'Xdata');
y=get(h1,'Ydata');
a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
if(number)
    text(h,x, y, c);
    h1.CData=[1 1 1];
end
m=max(max([x,y]));
axis([0 m 0 m])
refline(1)
xlabel('Bits Visual Cue Info')
ylabel('Bits Auditory Cue Info')
title('There is more visual info than auditory info in the population response')

h=subplot(2,2,3), scatter(Vinternal_choice_info,Ainternal_choice_info,ssize,cc(attend+1,:));
h1=get(h,'Children');
x=get(h1,'Xdata');
y=get(h1,'Ydata');
a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
if(number)
    text(h,x, y, c);
    h1.CData=[1 1 1];
end
m=max(max([x,y]));
axis([0 m 0 m])
refline(1)
xlabel('Bits Pop resp that drives chioce independent of Visual stim')
ylabel('Bits Pop resp that drives chioce independent of Auditory stim')
title('Visual information correlates more with behavior')

h=subplot(2,2,4), scatter(Vnon_readout_sensory_info,Anon_readout_sensory_info,ssize,cc(attend+1,:));
h1=get(h,'Children');
x=get(h1,'Xdata');
y=get(h1,'Ydata');
a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
if(number)
    text(h,x, y, c);
    h1.CData=[1 1 1];
end
m=max(max([x,y]));
axis([0 m 0 m])
refline(1)
xlabel('Non-Readout Visual Info')
ylabel('Non-Readout Auditory Info')
title('Alot of the Auditory information is not used to drive behavior')



% 
% MIjoint_r_s = - sum(mean(psgr).*log2(mean(psgr))) + mean(sum(psgr.*log2(psgr),2));
% MIjoint_r_b = - sum(mean(pbgr).*log2(mean(pbgr))) + mean(sum(pbgr.*log2(pbgr),2));
% 
% MI_r_sgb0= calcMI(alldata(resp==0,:),vstim(resp==0)>0);
% MI_r_sgb1= calcMI(alldata(resp==1,:),vstim(resp==1)>0);
% MI_r_sgb = mean(resp==0)*MI_r_sgb0 + mean(resp==1)*MI_r_sgb1;
% 
% MI_r_bgs0= calcMI(alldata(vstim==0,:),resp(vstim==0)>0);
% MI_r_bgs1= calcMI(alldata(vstim>0,:),resp(vstim>0)>0);
% MI_r_bgs = mean(vstim==0)*MI_r_bgs0 + mean(vstim>0)*MI_r_bgs1;
% 
% MI_s_b = calcMI(zscore(vstim>0),resp>0);
% MI_b_s = calcMI(zscore(resp>0),vstim>0);
% 
% MI_rs(n) = MI_r_s;
% MI_rb(n) = MI_r_b;
% MI_sb(n) = MI_s_b;
% MIjoint_rs(n) = MIjoint_r_s;
% MIjoint_rb(n) = MIjoint_r_b;
% MIjoint_sbr(n) = MI_r_sb;
% CoI_rsgb(n) = MIjoint_r_s - MI_r_sgb;
% CoI_rbgs(n) = MIjoint_r_b - MI_r_bgs;

