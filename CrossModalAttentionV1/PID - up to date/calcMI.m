function [MI,pc,conf,p,model] = calcMI(X,Y,sparse)
    if(~exist('sparse','var'))
        sparse=0;
    end
    if(islogical(Y))
        Y=Y+1;
    end
    X=X-mean(X);
    X=[X,ones(size(X,1),1)];
    model=VBMNLR(max(Y),size(X,2),sparse);
    model.update(Y,X,200);
    [Yhat,p]=model.getSimplePredictions(X);
    pY=mean(p);
    
    [m,loc]=max((p./pY)');  % Yhat is MAP estimage, this is the ML estimate
    pc=mean(loc'==Y);
%    pc=mean(Yhat==Y);

    idx=(p>0);
    pY=mean(p);
    MI = -sum(pY.*log2(pY)) + sum(sum(p(idx).*log2(p(idx))))/size(p,1);
    
    conf=zeros(size(p,2));
    [m,loc]=max((p./pY)');
    for i=1:size(p,2)
    for j=1:size(p,2)
        conf(i,j) = mean(loc(Y==j)==i);
    end
    end
        
end

