function [MI,pc,conf,p,stat] = calcMI(X,Y)
    if(islogical(Y))
        dv=mean(Y);
        parms=statset('glmfit');
        parms.MaxIter=200;
        B=glmfit(X,Y,'binomial');
        p = glmval(B,X,'logit');
        CD=sum(p(Y==1)>dv)/sum(Y==1);
        FP=sum(p(Y==0)>dv)/sum(Y==0);
        pc = (sum(p(Y==1)>dv) + sum(p(Y==0)<dv))/length(Y);
        pY=mean(Y);
        MI = - pY*log2(pY) - (1-pY)*log2(1-pY) ...
             + mean(p.*log2(p)) + mean((1-p).*log2(1-p));
        conf(1,1)=1-FP;conf(2,1)=FP;conf(1,2)=1-CD;conf(2,2)=CD;
    else
        B=mnrfit(X,Y);
        p=mnrval(B,X);
        pY=mean(p);
        conf=zeros(size(p,2));
        [m,loc]=max((p./pY)');
        for i=1:size(p,2)
        for j=1:size(p,2)
            conf(i,j) = mean(loc(Y==j)==i);
        end
        end
        pc=mean(loc'==Y);
        MI = -sum(pY.*log2(pY)) + mean(sum(p.*log2(p),2));
    end
    
end

