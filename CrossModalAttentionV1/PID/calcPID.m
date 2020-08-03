function[PID,MI,ENT,PR] = calcPID(D,pncs,binaryS,useclusters)
%
%Run outside of function
% addpath ./PID/
% addpath ./PID/glpkmex/
% addpath ./PID/glpkmex/win64/
%
%D{n} contains the data obtained from preprocess.m
%
%
% pncs is the number of clusters to cluster p(s,b|r) for use with pca data and mnlr fit.
%
% binaryS is a logical which if true (or 1) converts the stimulus to 1 for
% all values greater than 1.
%
% useclusters is a logical which if true (or 1) uses D{n}.
% PID contains the information theoretic quantites
%
% MI an independent computation of pairwise mutual information
% ENT is the entropy
%
for n=1:length(D)
if(binaryS)
    D{n}.Vis.S=D{n}.Vis.S>0;
    D{n}.Aud.S=D{n}.Aud.S>0;
end
%
%Directly compute mutual information between variabel pairs.    
if(~useclusters)
    MI.Vis_RS(n) = calcMI(D{n}.Vis.R,D{n}.Vis.S+1);
    MI.Vis_SB(n) = calcMI(D{n}.Vis.S,D{n}.Vis.B+1);
    MI.Vis_RB(n) = calcMI(D{n}.Vis.R,D{n}.Vis.B+1);

    MI.Aud_RS(n) = calcMI(D{n}.Aud.R,D{n}.Aud.S+1);
    MI.Aud_SB(n) = calcMI(D{n}.Aud.S,D{n}.Aud.B+1);
    MI.Aud_RB(n) = calcMI(D{n}.Aud.R,D{n}.Aud.B+1);
else
    MI=NaN;
end

[C,IA,IC]=unique(D{n}.Vis.S);
pIA=IA/sum(IA);
ENT.Vis_S(n)=-sum(pIA.*log2(pIA));

[C,IA,IC]=unique(D{n}.Aud.S);
pIA=IA/sum(IA);
ENT.Aud_S(n)=-sum(pIA.*log2(pIA));

[C,IA,IC]=unique(D{n}.Vis.B);
pIA=IA/sum(IA);
ENT.Vis_B(n)=-sum(pIA.*log2(pIA));

[C,IA,IC]=unique(D{n}.Aud.B);
pIA=IA/sum(IA);
ENT.Aud_B(n)=-sum(pIA.*log2(pIA));

if(~useclusters)
    ENT.Aud_R(n)=NaN;
    ENT.Vis_R(n)=NaN;
    
    ns=length(unique(D{n}.Vis.S));
    [MI.Vis.RSnB(n),~,~,pvbgr] = calcMI(D{n}.Vis.R,D{n}.Vis.S+ns*D{n}.Vis.B+1);
    m=max(D{n}.Vis.S+ns*D{n}.Vis.B+1);
    if(m<2*ns)
        pvbgr(:,m+1:2*ns)=0;
    end

    ns=length(unique(D{n}.Aud.S));
    [MI.Aud_RSnB(n),~,~,pabgr] = calcMI(D{n}.Aud.R,D{n}.Aud.S+ns*D{n}.Aud.B+1);
    m=max(D{n}.Aud.S+ns*D{n}.Aud.B+1);
    if(m<2*ns)
        pabgr(:,m+1:2*ns)=0;
    end

    [z,c]=kmeans(pvbgr,pncs);
    zc=hist(z,pncs);
    zc=zc/sum(zc);
    ns=length(unique(D{n}.Vis.S));
    pv=zeros(ns,pncs,2);
    for s=1:ns
    for b=1:2
    for r=1:pncs    
        pv(s,r,b)=sum(pvbgr(z==r,s+ns*(b-1)));
    end
    end
    end

    [z,c]=kmeans(pabgr,pncs);
    zcs=hist(z,pncs);
    zc=zc/sum(zc);
    ns=length(unique(D{n}.Aud.S));
    pa=zeros(ns,pncs,2);
    for s=1:ns
    for b=1:2
    for r=1:pncs
        pa(s,r,b)=sum(pabgr(z==r,s+ns*(b-1)));
    end
    end
    end
else %use clusters
    ns=length(unique(D{n}.Vis.S));
%     r1=unique(D{n}.Vis.Rlabels);
%     nr=length(r1);
    nr=size(D{n}.clusters,1);
    pv=zeros(ns,nr,2);
    for s=1:ns
    for r=1:nr
    for b=1:2
%        pv(s,r,b)=sum(D{n}.Vis.S+1==s & D{n}.Vis.Rlabels==r1(r) & D{n}.Vis.B+1==b);
        pv(s,r,b)=sum(D{n}.Vis.S+1==s & D{n}.Vis.Rlabels==r & D{n}.Vis.B+1==b);
    end
    end
    end
    
    ns=length(unique(D{n}.Aud.S));
%    rl=unique(D{n}.Aud.Rlabels);
%    nr=length(r1);
    nr=size(D{n}.clusters,1);
    pa=zeros(ns,nr,2);
    for s=1:ns
    for r=1:nr
    for b=1:2
%        pa(s,r,b)=sum(D{n}.Aud.S+1==s & D{n}.Aud.Rlabels==rl(r) & D{n}.Aud.B+1==b);
        pa(s,r,b)=sum(D{n}.Aud.S+1==s & D{n}.Aud.Rlabels==r & D{n}.Aud.B+1==b);
    end
    end
    end
        
    [C,IA,IC]=unique(D{n}.Aud.Rlabels);
    pIA=IA/sum(IA);
    ENT.Aud_R(n)=-sum(pIA.*log2(pIA));
    
    [C,IA,IC]=unique(D{n}.Vis.Rlabels);
    pIA=IA/sum(IA);
    ENT.Vis_R(n)=-sum(pIA.*log2(pIA));
    
end
pa=pa/sum(pa(:));
pv=pv/sum(pv(:));

PR.sbr.Aud{n}=permute(pa,[1 3 2]);
PR.sbr.Vis{n}=permute(pv,[1 3 2]);

stemp=sum(sum(PR.sbr.Aud{n},2),1);
stemp(stemp==0)=1;
PR.sbgr.Aud{n}=PR.sbr.Aud{n}./stemp;

stemp=sum(sum(PR.sbr.Vis{n},2),1);
stemp(stemp==0)=1;
PR.sbgr.Vis{n}=PR.sbr.Vis{n}./stemp;

PR.sgr.Aud{n}=squeeze(sum(PR.sbgr.Aud{n},2));
PR.sgr.Vis{n}=squeeze(sum(PR.sbgr.Vis{n},2));
PR.bgr.Aud{n}=squeeze(sum(PR.sbgr.Aud{n},1));
PR.bgr.Vis{n}=squeeze(sum(PR.sbgr.Vis{n},1));

PR.r.Vis{n}=squeeze(sum(sum(pv,3),1));
PR.r.Aud{n}=squeeze(sum(sum(pa,3),1));
tic
[PID.II.Vis(n), PID.MI_SR.Vis(n), PID.MI_BR.Vis(n), PID.MI_SB.Vis(n), PID.NRSI.Vis(n), PID.ICI.Vis(n), PID.MISB_not_in_R(n).Vis] ...
=intersection_information(pv);

[PID.II.Aud(n), PID.MI_SR.Aud(n), PID.MI_BR.Aud(n), PID.MI_SB.Aud(n), PID.NRSI.Aud(n), PID.ICI.Aud(n), PID.MISB_not_in_R(n).Aud] ...
=intersection_information(pa);
[n/length(D), n, toc]

end

end

