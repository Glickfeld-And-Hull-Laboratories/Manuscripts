
function[PID,MI,ENT] = calcPID(D,npc,binaryS)

%Run outside of function
%addpath ./PID/
%addpath ./PID/glpkmex/
%addpath ./PID/glpkmex/win64/

%D{n} contains the data obtained from preprocess.m
%
% npc is the number of clusters to put posterior probabilities into.  This
% is needed to speed up the PID calculation  A good number seems to be
% about 100.
%
% PID contains the information theoretic quantites


% MI an independent computation of pairwise mutual information
% ENT is the entropy

for n=1:length(D)
if(binaryS)
    D{n}.Vis.S=double(D{n}.Vis.S>0);
    D{n}.Aud.S=double(D{n}.Aud.S>0);
end

%Directly compute mutual information between variabel pairs.    
MI.Vis_RS(n) = calcMI(D{n}.Vis.R,D{n}.Vis.S+1);
MI.Vis_SB(n) = calcMI(D{n}.Vis.S,D{n}.Vis.B+1);
MI.Vis_RB(n) = calcMI(D{n}.Vis.R,D{n}.Vis.B+1);

MI.Aud_RS(n) = calcMI(D{n}.Aud.R,D{n}.Aud.S+1);
MI.Aud_SB(n) = calcMI(D{n}.Aud.S,D{n}.Aud.B+1);
MI.Aud_RB(n) = calcMI(D{n}.Aud.R,D{n}.Aud.B+1);



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

[z,c]=kmeans(pvbgr,npc);
zc=hist(z,npc);
zc=zc/sum(zc);
ns=length(unique(D{n}.Vis.S));
pv=zeros(ns,npc,2);
for s=1:ns
for b=1:2
for r=1:npc    
    pv(s,r,b)=sum(pvbgr(z==r,s+ns*(b-1)));
end
end
end

[z,c]=kmeans(pabgr,npc);
zc=hist(z,npc);
zc=zc/sum(zc);
ns=length(unique(D{n}.Aud.S));
pa=zeros(ns,npc,2);
for s=1:ns
for b=1:2
for r=1:npc
    pa(s,r,b)=sum(pabgr(z==r,s+ns*(b-1)));
end
end
end


[PID.II.Vis(n), PID.MI_SR.Vis(n), PID.MI_BR.Vis(n), PID.MI_SB.Vis(n), PID.NRSI.Vis(n), PID.ICI.Vis(n), PID.MISB_not_in_R(n).Vis] ...
=intersection_information(pv/sum(pv(:)));

[PID.II.Aud(n), PID.MI_SR.Aud(n), PID.MI_BR.Aud(n), PID.MI_SB.Aud(n), PID.NRSI.Aud(n), PID.ICI.Aud(n), PID.MISB_not_in_R(n).Aud] ...
=intersection_information(pa/sum(pa(:)));

end

end

