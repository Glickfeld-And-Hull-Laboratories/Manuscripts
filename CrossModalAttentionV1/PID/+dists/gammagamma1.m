classdef gammagamma1 < handle
    properties
        % Gamma(alpha,alpha)
        % alpha hyperparameters
        a_0
        b_0
        
        % alpha model params
        a
        b

        % variational parameter
        xi
        
        % assumes that x is Gamma(alpha,alpha) and utilizes a linear
        % lower bound to the convex function f = x*log(x) - gammaln(x) - log(x);
    end
    
    methods
        function self = gammagamma1(a_0,b_0)
            self.a_0 = a_0;
            self.b_0 = b_0;
            self.a = a_0;%.*(1+rand(size(a_0)));
            self.b = b_0;%.*(1+rand(size(a_0)));

            self.xi = self.a./self.b;
            
        end
        
        function update(self,Ex,Elogx,n)  
            for i=1:5  % run twice so that the speed of the update is comparable to 
                       % updates for conjugate priors.
                self.xi = self.alphamean;
                self.a = self.a_0 + n;
                self.b = self.b_0 + n.*Ex - n.*self.fprime  - Elogx.*n;
%                self.b = min(self.b,self.b_0);
                self.xi = self.alphamean;
            end
        end
        
        function updateSS(self,Ex,Elogx,n)  
            if(n>0)
                for i=1:5  % run twice so that the speed of the update is comparable to 
                           % updates for conjugate priors.
                    self.xi = self.alphamean;
                    self.a = self.a_0 + n;
                    self.b = self.b_0 + Ex.*n - n.*self.fprime - Elogx.*n;
                    self.xi = self.alphamean;
                end
            else
                self.a=self.a_0;
                self.b=self.b_0;
                self.xi = self.alphamean;
                
            end
        end
        
        function rawupdate(self,data,p)
            if(~exist('p','var'))
               p=ones(size(data,1),1);
            end
            idx=find(~isnan(sum(data,2)));
            n=sum(p(idx));
            SEx = p(idx)'*data(idx,:);
            SElogx = p(idx)'*log(data(idx,:));
            self.updateSS(SEx',SElogx',n);
        end
        
        function res = mean(self)
            res=ones(size(self.a));
        end
        
        function res = meaninv(self)
            res=self.a./(self.a-1);
        end
        
        function res = var(self)
            if(self.a>1)
                res=self.b./(self.a-1);
            else
                res = Inf;
            end        
        end

        function res = secondmoment(self)
            if(self.a>1)
                res = self.b./(self.a-1) + 1;
            else
                res = Inf;
            end        
        end

        function res = alphamean(self)
            res = self.a ./ self.b;
        end
        
        function res = alphaloggeomean(self)
            res = psi(self.a) - log(self.b);
        end
                        
        function res = Ealpha(self)
            res = self.a ./ self.b;
        end
        
        function res = Elogalpha(self)
            res = psi(self.a) - log(self.b);
        end
                
        function res = KLqprior(self) 
            res = (self.a-self.a_0).*psi(self.a) - gammaln(self.a) + gammaln(self.a_0) ...
                + self.a_0.*(log(self.b)-log(self.b_0)) + self.a.*(self.b_0./self.b-1);
        end

        function res = lowerboundcontrib(self)
            res = - self.KLqprior;
        end
        
%         function res = Eloglikelihood(self,data) % assumes gamma likelihood
%             res = log(data)*(self.alphamean-1) - data*(self.alphamean) ...
%                 + sum(self.alphamean.*self.betaloggeomean) ... 
%                 + sum(self.alphaloggeomean) ...
%                 + sum(self.f) + sum(self.fprime.*(self.alphamean-self.xi));
%             res(isnan(res)) = 0;
%         end
        
%         function res = EloglikelihoodSS(self,Ex,Elogx,n) % NOT VALIDATED
%             res = Elogx*(self.alphamean-1) - Ex*(self.alphamean.*self.betamean) ...
%                 + sum(self.alphamean.*self.betaloggeomean) ... 
%                 + sum(self.alphaloggeomean) ...
%                 + sum(self.f) + sum(self.fprime.*(self.alphamean-self.xi));
%             res=res*n;
%             res(isnan(res)) = 0;
%         end
        
        function res = logZ(self)
            res = self.f + self.Elogalpha;
        end
        
        function res = logZp(self)
            xisave=self.xi;
            asave=self.a;
            bsave=self.b;
            
            self.xi=self.a_0./self.b_0;
            self.a=self.a_0;
            self.b=self.b_0;
            
            res = self.logZ;
            
            self.xi=xisave;
            self.a=asave;
            self.b=bsave;
        end
        
        function L = fit(self,data,tol,maxiters)
            Ex=mean(data);
            Elogx=mean(log(data));
            [N,dim] = size(data);
            
            self.a_0 = ones(1,dim);
            self.b_0 = ones(1,dim);
            self.a = (1+rand(1,dim));
            self.b = (1+rand(1,dim));
            self.xi = self.a_0./self.b_0;
            k=0;
            Llast=-Inf;
            self.update(Ex,Elogx,N);

            L=sum(self.Eloglikelihood(data)) - sum(self.KLqprior);
            while( L-Llast>tol*abs(L) && k<maxiters )
                k=k+1;
                self.update(Ex,Elogx,N);
                Llast=L;
                L=sum(self.Eloglikelihood(data)) - sum(self.KLqprior); 
            end
            
            if(k==maxiters)
                'maxiters reached'
            else
               ['converged in ',num2str(k),' iterations to an ELBO of ',num2str(L)] 
            end
        end
        
        function res = f(self)
            res = self.xi.*log(self.xi) - gammaln(self.xi) - log(self.xi);
        end
               
        function res = fprime(self)
            res = log(self.xi) + 1 - psi(self.xi) - 1./self.xi;
        end
        
    end
end




