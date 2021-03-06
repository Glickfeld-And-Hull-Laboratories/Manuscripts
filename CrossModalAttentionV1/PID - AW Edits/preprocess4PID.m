%%% ATTENDED VS UNATTENDED 
%%% trial type cued by presence or absense of tone.  
%%% invalid visual trial means auditory target comes
%%% invalid auditory trial means visual comes first
%%% 
%%% a mouse has attention if an invalid stimulus does not elilct a response
%%% (as often).    These were rare.  

% clear all
% close all
% load('attnData_neuronResp_trialOutcome.mat')
% 
% npcs=15;  %maxumum number of principle components to save
% nc=60;    %number of clusters of patterns (kmeans)
function D = preprocess4PID(data,npcs,nc,binaryS)
%%% set npcs equal to 'expt' to calculate the number of pcs to be 10% of
%%% the total trials (max 15 pcs)




%Add frame number as a categorical regressor for predicting behavior.

for n=1:length(data)

attend(n)=data(n).hasAttention;
NT=length(data(n).dff);

% option to determine number of PCs used 
    if strcmp(npcs,'expt')
        NPC = round(0.1*NT);
    elseif strcmp(npcs,'expt_plus')
        NPC = round(0.1*NT);
        if NPC > 20 && binaryS
            NPC = 20;
        elseif NPC > 15 && ~binaryS
            NPC = 15;
        end
    elseif isnumeric(npcs)
        NPC = npcs;
    else
        error("set npcs to number or expt or expt_plus")
    end

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

%%%% ASHLEY LOOK HERE!!!!
% vstim(vstim>0)=1;
% astim(astim>0)=1;
vstim(vstim>4)=4;
idx = frame>3 & cue==1;% & vstim<=5; 


cue1=cue(idx);
vstim1=vstim(idx); 
astim1=astim(idx);
frame1=frame(idx);
resp1=resp(idx);
alldata1=alldata(idx,:);

%%%%%% ASHLEY LOOK HERE TOO (ITS VERY SIMILIAR)
astim(astim>4)=4;
idx = frame>3 & cue==0 ;%& astim<=4; %& cue==1;% & vstim<=3;

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
idx1=[ones(size(alldata1,1),1);zeros(size(alldata2,1),1);]==1;
idx2=~idx1;

%Cluster reduction
[coeff,score]=pca(alldata);
score=score(:,1:min(NPC,size(alldata,2)));
[z,c]=kmeans(score,nc);
datatemp=zeros(size(alldata,1),nc);
for i=1:size(alldata,1)
    datatemp(i,z(i))=1;
end
idx=sum(datatemp,1)>0;
datatemp=datatemp(:,idx);
alldata=datatemp;
alldata1=alldata(idx1,:);
alldata2=alldata(idx2,:);
idx=sum(alldata1,1)>1;
alldata1=alldata1(:,idx);
idx=sum(alldata2,1)>1;
alldata2=alldata2(:,idx);

D{n}.clusters=c;
D{n}.Rlabels=z;
D{n}.Vis.Rlabels=z(idx1);
D{n}.Aud.Rlabels=z(idx2);
% 
%%%PCA reduction

mu=mean(alldata);
[coeff,score,lambda]=pca(alldata);
score=score(:,1:min(NPC,size(alldata,2)));

alldata=zscore(score);

D{n}.mu=mu;
D{n}.coeff=coeff;
D{n}.lambda=lambda;

D{n}.attend=data(n).hasAttention;
D{n}.Vis.Rraw=alldata1;
D{n}.Aud.Rraw=alldata2;

alldata1=alldata(idx1,:);
alldata2=alldata(idx2,:);

D{n}.Vis.B=resp1;
D{n}.Vis.R=alldata1;
D{n}.Vis.S=vstim1;
D{n}.Vis.NS=astim1;
D{n}.Vis.frame=frame1;

D{n}.Aud.B=resp2;
D{n}.Aud.R=alldata2;
D{n}.Aud.S=astim2;
D{n}.Aud.NS=vstim2;
D{n}.Aud.frame=frame2;


end
%%%end proprocess
end