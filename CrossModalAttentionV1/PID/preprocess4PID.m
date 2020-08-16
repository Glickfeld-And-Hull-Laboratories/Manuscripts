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
% addpath ./PID/
% addpath ./PID/glpkmex/
% addpath ./PID/glpkmex/win64/
function D = preprocess4PID(data,npcs,nc,binaryS)
% npcs=20;  %maxumum number of principle components to save
% nc=80;    %number of clusters of patterns (kmeans)
mergstims=2;
startframe=3;
nopca4clusters=1;
zscorepcs=1;
for n=1:length(data);
    % option to determine number of PCs used 
    NT=length(data(n).dff);
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
%%%% insert preprocess loop here
%Add frame number as a categorical regressor for predicting behavior.
attend(n)=data(n).hasAttention;
NT=length(data(n).dff);

alldata=[];
resp=[];
vstim=[];
astim=[];
frame=[];
cue=[];
islast2frames=[];
for t=1:NT
    if(sum(data(n).dff{t},2)==0)
    else
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
        if(np>=2)
            islast2frames(length(frame)-1,1)=1;
        end
        islast2frames(length(frame),1)=1;
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
        if(np>=2)
            islast2frames(length(frame)-1,1)=1;
        end
        islast2frames(length(frame),1)=1;
        cue=[cue;zeros(np,1);];%0 for auditory trials

    end
    end    
end
vstim(isnan(vstim))=0;
[C,IA,IC]=unique(vstim);
vstim=IC-1;  % task difficulty 0 = distractor and 1 to 6 is 16 to 90 degrees.


astim(isnan(astim))=0;
[C,IA,IC]=unique(astim);
astim=IC-1; % no sound == 0;

%%%% ASHLEY LOOK HERE!!!!

% vstim(vstim>0)=1;
% astim(astim>0)=1;
%vstim(vstim>4)=4;
if(mergstims==1)
    vstim(vstim>4)=4;
    astim(astim>4)=4;
elseif(mergstims==2)
    vstim(vstim==2|vstim==3)=2;
    vstim(vstim==4|vstim==5)=3;
    vstim(vstim>5)=4;
    astim(astim==2|astim==3)=2;
    astim(astim==4|astim==5)=3;
    astim(astim>5)=4;
elseif(mergstims==3)
    vstim(vstim==1|vstim==2)=1;
    vstim(vstim==3|vstim==4)=2;
    vstim(vstim>4)=3;
    astim(astim==1|astim==2)=1;
    astim(astim==3|astim==4)=2;
    astim(astim>4)=3;    
elseif(mergstims==4)
    vstim(vstim==1|vstim==2)=1;
    vstim(vstim==3|vstim==4)=2;
    vstim(vstim==5|vstim==6)=3;
    vstim(vstim>6)=4;    
    astim(astim==1|astim==2)=1;
    astim(astim==3|astim==4)=2;
    astim(astim==5|astim==6)=3;
    astim(astim>6)=4;
elseif(mergstims==5)
    vstim(vstim==4|vstim==5)=4;
    vstim(vstim>=6)=5;    
    astim(astim==4|astim==5)=4;
    astim(astim>=6)=5;
end
idx = frame>=startframe & cue==1 & vstim<=4 & islast2frames==1; 

cue1=cue(idx);
vstim1=vstim(idx); 
astim1=astim(idx);
frame1=frame(idx);
resp1=resp(idx);
alldata1=alldata(idx,:);
islast2frames1=islast2frames(idx);
%%%%%% ASHLEY LOOK HERE TOO (ITS VERY SIMILIAR)
%astim(astim>4)=4;

idx = frame>=startframe & cue==0 & astim<=4 & islast2frames==1;
cue2=cue(idx);
vstim2=vstim(idx); 
astim2=astim(idx);
frame2=frame(idx)-1;
resp2=resp(idx);
alldata2=alldata(idx,:);
islast2frames2=islast2frames(idx);

cue=[cue1;cue2];
vstim=[vstim1;vstim2];
astim=[astim1;astim2];
frame=[frame1;frame2];
islast2frames=[islast2frames1;islast2frames2];
resp=[resp1;resp2];
alldata=[alldata1;alldata2;];
idx1=[ones(size(alldata1,1),1);zeros(size(alldata2,1),1);]==1;
idx2=~idx1;

%PCA reduction

[coeff,score,lambda]=pca(alldata);
score=score(:,1:min(npcs,size(alldata,2)));
mu=mean(alldata);
alldata1=alldata(idx1,:);
alldata2=alldata(idx2,:);
D{n}.Vis.Rraw=alldata1;
D{n}.Aud.Rraw=alldata2;
D{n}.mu=mu;
D{n}.coeff=coeff;
D{n}.lambda=lambda;

if(nopca4clusters)
    [z,c]=kmeans(alldata,nc);
else
    [z,c]=kmeans(score,nc);
end
D{n}.clusters=c;
D{n}.Rlabels=z;
D{n}.Vis.Rlabels=z(idx1);
D{n}.Aud.Rlabels=z(idx2);

if(zscorepcs)
    score=zscore(score);
end

D{n}.attend=data(n).hasAttention;
D{n}.Vis.B=resp1;
D{n}.Vis.R=score(idx1,:);
D{n}.Vis.S=vstim1;
D{n}.Vis.NS=astim1;
D{n}.Vis.frame=frame1;
D{n}.Vis.islast2frames=islast2frames1;

D{n}.Aud.B=resp2;
D{n}.Aud.R=score(idx2,:);
D{n}.Aud.S=astim2;
D{n}.Aud.NS=vstim2;
D{n}.Aud.frame=frame2;
D{n}.Aud.islast2frames=islast2frames2;

end
%%%end proprocess
end