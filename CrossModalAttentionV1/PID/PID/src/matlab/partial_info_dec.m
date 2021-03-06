function [I_shar, I_syn, I_unx, I_uny, q_opt] = partial_info_dec(p)
    % partial_info_dec(p) Partial Information Decomposition of a trivariate probability distribution.
    % 
    % [I_shar, I_syn, I_unx, I_uny] = partial_info_dec(p) gives the four atoms of the
    % partial information decomposition following the proposal by
    % Bertschinger et al. (2014).
    %
    % p must be a three-dimensional array representing the joint
    % probability distribution p(x,y,z) of the three discrete
    % variables X, Y and Z. Z is the target variable. The values
    % assumed by X, Y and Z are implicitly mapped to indices such
    % that, for instance, p(1,3,2) is the probability of X=1,Y=3,Z=2,
    % and so on.
    %
    % [I_shar, I_syn, I_unx, I_uny] are the atoms of the PID:
    %
    % I_shar = SI(Z:{X;Y})  (shared information, or redundancy, between X
    % and Y about X)
    %
    % I_syn = CI(Z:{X,Y})   (complementary information, or synergy, between
    % X and Y about Z)
    % 
    % I_unx = UI(Z:{X\Y})   (unique information in X about Z)
    %
    % I_uny = UI(Z:{Y\Z})   (unique information in Y about Z)
    %
    % q_opt is the probability distribution that solves the optimization
    % problem defined in Bertschinger et al 2014, which is used to compute
    % the quantities above.
    
    [dimx, dimy, dimz] = size(p);

    % hardcoded parameters - only change if you know what you're doing
    accuracy = 0.001;
    glpk_verbosity = 0; % set to 0 to suppress all output from glpk.
    line_search=0; % the increment is optimized with a line-search
                   % only if necessary: the code below
                   % automatically sets line_search=1 when needed
    iter_limit = 10^4;
    
    if dimx==1
        
        p23 = reshape(p(1,:,:),dimy,dimz);
        I_yz = p23 .* log2(p23 ./ repmat(sum(p23), [dimy 1]) ./ repmat(sum(p23,2), [1 dimz]));
        I_yz = sum(I_yz(p23 > 0));
        
        I_shar=0;
        I_syn=0;
        I_unx=0;
        I_uny=I_yz;
        
        q_opt = p;
        
    elseif dimy==1
        
        p13 = reshape(p(:,1,:),dimx,dimz);
        I_xz = p13 .* log2(p13 ./ repmat(sum(p13), [dimx 1]) ./ repmat(sum(p13,2), [1 dimz]));
        I_xz = sum(I_xz(p13 > 0));
        
        I_shar=0;
        I_syn=0;
        I_unx=I_xz;
        I_uny=0;
        
        q_opt = p;
        
    else
        
        GAMMA=zeros(dimx,dimy,dimz,dimx-1,dimy-1,dimz); % this is the concatenation of all the (dimx-1)*(dimy-1)*(dimz) Gamma matrices defined in Bertschinger2014, each of which has
        %the same dimensions (dimx, dimy, dimz) as the input p.
        
        for zz=1:dimz
            for xx=1:dimx-1
                for yy=1:dimy-1
                    GAMMA(xx,yy,zz,xx,yy,zz) = 1;
                    GAMMA(xx+1,yy,zz,xx,yy,zz) = -1;
                    GAMMA(xx,yy+1,zz,xx,yy,zz) = -1;
                    GAMMA(xx+1,yy+1,zz,xx,yy,zz) = 1;
                end
            end
        end
        
        
        check=0;% when check==1, the algorithm is over and the output is returned
        
        parameters=zeros(1,ceil(dimz*(dimx-1)*(dimy-1))); % the coefficients of the matrix q in the basis of the Gamma matrices
        
        iter=0; % counts the iterations of the algorithm
        
        q=p; % the starting point of the algorithm is trivially set to the input p. Different, smarter starting points seem to 
        % make little difference in terms of performance.
        
        coeff_prev=parameters';% when iter=0, the coefficients of the iteration -1 are trivially set equal to zero too.
        
        while check==0 % iteration loop
            
            q(q<0)=0;% eliminate tiny negative entries in q which result from the limited numerical precision
            
            % I_q(X:Y)
            q12 = sum(q, 3);
            I_xy = q12 .* log2(q12 ./ repmat(sum(q12), [dimx 1]) ./ repmat(sum(q12,2), [1 dimy]));
            I_xy = sum(I_xy(q12 > 0));
            
            % I_q(X:Y|Z)
            I_cond_xy_z = q .* log2(q ./ repmat(sum(sum(q), 2), [dimx dimy 1]) ./ ...
                ( repmat(sum(q,2), [1 dimy 1]) ./ repmat(sum(sum(q), 2), [dimx dimy 1]) .* ...
                repmat(sum(q),   [dimx 1 1]) ./ repmat(sum(sum(q), 2), [dimx dimy 1]) ) );
            I_cond_xy_z = sum(I_cond_xy_z(q > 0));
            
            % update coI_q
            co_I = I_xy - I_cond_xy_z;
            
            if iter==0
                co_I_in=co_I;
            end
            
            %Franke-Wolf optimization algorithm
            %1) determine search direction
            
            %calculating the gradient: we have an analytical expression of
            %the gradient of the object function
            
            deriv = log2(q(1:dimx-1,1:dimy-1,:) .* q(2:dimx,2:dimy,:)) - ...
                log2(q(1:dimx-1,2:dimy,:) .* q(2:dimx,1:dimy-1,:)) + ...
                log2( repmat(sum(q(1:dimx-1,2:dimy,:), 3), [1 1 dimz]) .* repmat(sum(q(2:dimx,1:dimy-1,:), 3), [1 1 dimz]) ) - ...
                log2( repmat(sum(q(1:dimx-1,1:dimy-1,:), 3), [1 1 dimz]) .* repmat(sum(q(2:dimx,2:dimy,:), 3), [1 1 dimz]) );
            
            %get rid of nonsense values of deriv coming from finite
            %numerical precision
            deriv(isnan(deriv))=0;
            deriv(isinf(deriv))=0;
            
            %store the coefficients of the q of the current iteration. For
            %each value of z there is a coeff, then coeff_tot concatenates
            %all coeffs.
            coeff_tot=zeros(ceil((dimx-1)*(dimy-1)),dimz);
            
            %the optimization can be formally divided into dimz
            %optimizations, one for each value of the variable Z. This loop
            %could thus be parallelized
            for zz=1:dimz
                
                %the constraints on q, for each z, are implemented via the inequality A*coeff<=b
                b=zeros(ceil(2*dimx*dimy),1);
                A=zeros(ceil(2*dimx*dimy),ceil((dimx-1)*(dimy-1)));
                
                %the number of constraints to be imposed
                count=1;
                
                for xx=1:dimx-1
                    for yy=1:dimy-1
                        
                        % set the constraints for the top-left entries q(1,1,zz)
                        if (xx==1 && yy==1)
                            
                            A(count,(floor(xx)-1)*(dimy-1)+yy)=-1;
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy)=1;
                            
                            b(count)=p(xx,yy,zz);
                            b(count+ceil(dimx*dimy))=1-p(xx,yy,zz);
                            
                            count=count+1;
                        end
                        
                        % bottom right entries
                        if (xx==dimx-1 && yy==dimy-1)
                            
                            A(count,(floor(+xx)-1)*(dimy-1)+yy)=-1;
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy)=1;
                            
                            b(count)=p(xx+1,yy+1,zz);
                            b(count+ceil(dimx*dimy))=1-p(xx+1,yy+1,zz);
                            
                            count=count+1;
                        end
                        
                        % top right entries
                        if (xx==1 && yy==dimy-1)
                            
                            A(count,(floor(+xx)-1)*(dimy-1)+yy)=1;
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy)=-1;
                            
                            b(count)=p(xx,yy+1,zz);
                            b(count+ceil(dimx*dimy))=1-p(xx,yy+1,zz);
                            
                            count=count+1;
                        end
                        
                        %bottom left entries
                        
                        if (xx==dimx-1 && yy==1)
                            
                            A(count,(floor(+xx)-1)*(dimy-1)+yy)=1;
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy)=-1;
                            
                            b(count)=p(xx+1,yy,zz);
                            b(count+ceil(dimx*dimy))=1-p(xx+1,yy,zz);
                            
                            count=count+1;
                        end
                        
                        % top row entries, from left to right
                        if xx==1 && dimy>2 && yy<dimy-1
                            
                            A(count,(floor(+xx)-1)*(dimy-1)+yy)=1;
                            A(count,(floor(+xx)-1)*(dimy-1)+yy+1)=-1;
                            
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy)=-1;
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy+1)=1;
                            
                            
                            b(count)=p(xx,yy+1,zz);
                            b(count+ceil(dimx*dimy))=1-p(xx,yy+1,zz);
                            count=count+1;
                            
                        end
                        
                        % bottom row entries, from right to left
                        if xx==dimx-1 && dimy>2 && yy<dimy-1
                            
                            A(count,(floor(+xx)-1)*(dimy-1)+yy)=-1;
                            A(count,(floor(+xx)-1)*(dimy-1)+yy+1)=1;
                            
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy)=1;
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy+1)=-1;
                            
                            b(count)=p(xx+1,yy+1,zz);
                            b(count+ceil(dimx*dimy))=1-p(xx+1,yy+1,zz);
                            count=count+1;
                        end
                        
                        % left-most column, top-down
                        if yy==1 && dimx>2 && xx<dimx-1
                            
                            A(count,(floor(+xx)-1)*(dimy-1)+yy)=1;
                            A(count,(floor(+xx+1)-1)*(dimy-1)+yy)=-1;
                            
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy)=-1;
                            A(count+ceil(dimx*dimy),(floor(+xx+1)-1)*(dimy-1)+yy)=1;
                            
                            b(count)=p(xx+1,yy,zz);
                            b(count+ceil(dimx*dimy))=1-p(xx+1,yy,zz);
                            count=count+1;
                            
                        end
                        
                        % right-most column, top-down
                        if yy==dimy-1 && dimx>2 && xx<dimx-1
                            
                            A(count,(floor(+xx)-1)*(dimy-1)+yy)=-1;
                            A(count,(floor(+xx+1)-1)*(dimy-1)+yy)=1;
                            
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy)=1;
                            A(count+ceil(dimx*dimy),(floor(+xx+1)-1)*(dimy-1)+yy)=-1;
                            
                            b(count)=p(xx+1,yy+1,zz);
                            b(count+ceil(dimx*dimy))=1-p(xx+1,yy+1,zz);
                            count=count+1;
                        end
                        
                        % internal
                        if dimx>2 && dimy>2 && xx>1 && yy>1 && xx<=dimx-1 && yy<=dimy-1
                            
                            A(count,(floor(+xx)-1)*(dimy-1)+yy)=-1;
                            A(count,(floor(+xx-1)-1)*(dimy-1)+yy)=1;
                            A(count,(floor(+xx)-1)*(dimy-1)+yy-1)=1;
                            A(count,(floor(+xx-1)-1)*(dimy-1)+yy-1)=-1;
                            
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy)=1;
                            A(count+ceil(dimx*dimy),(floor(+xx-1)-1)*(dimy-1)+yy)=-1;
                            A(count+ceil(dimx*dimy),(floor(+xx)-1)*(dimy-1)+yy-1)=-1;
                            A(count+ceil(dimx*dimy),(floor(+xx-1)-1)*(dimy-1)+yy-1)=1;
                            
                            b(count)=p(xx,yy,zz);
                            b(count+ceil(dimx*dimy))=1-p(xx,yy,zz);
                            count=count+1;
                            
                        end
                        
                    end
                end
                
                %extract the portion of deriv pertaining to zz and adapt deriv_zz to the multiplication deriv_zz*coeff
                deriv_zz=permute(deriv(:,:,zz),[2 1 3]);
                deriv_zz=deriv_zz(:)';
                
                param = struct();
                param.lpsolver=2; % Interior point method. The simplex method called with param.lpsolver=1 or 3 doesn't work as well.
                param.msglev=glpk_verbosity;
                param.itcnt=10^10;
                param.mpsobj=0;
                param.mpsinfo=0;
                
                coeff = glpk(deriv_zz, A, b, [], [], repmat('U',1,length(b)), repmat('C',1,length(deriv_zz)), 1,param);
                
                coeff_tot(:,zz)=coeff;% update coeff_tot with the coeff from the current zz
                
            end
            
            coeff_tot=coeff_tot(:);
            
            p_k=p;
            
            for ind_gamma=1:ceil(dimz*(dimx-1)*(dimy-1))
                
                xx=mod(ceil(ind_gamma/(dimy-1))-1,dimx-1)+1;
                yy=mod(ind_gamma-1,dimy-1)+1;
                zz=ceil(ind_gamma/((dimx-1)*(dimy-1)));
                
                p_k=p_k+coeff_tot(ind_gamma,1).*GAMMA(:,:,:,xx,yy,zz);
                
            end
            
            deriv=permute(deriv,[2 1 3]);
            deriv=deriv(:)';
            
            %set the stopping criterion based on the duality gap, see Stratos;
            %iter must be larger than 1 because sometimes deriv takes 2 iters to
            %become different from zero, which is always its initial value.
            if (iter>1 && (dot(deriv,coeff_prev-coeff_tot)<=accuracy)) ...
                    || iter>iter_limit
                
                check=1; %exit the algorithm
                q_opt=q; % output the optimal distribution
                
                if co_I<co_I_in-accuracy
                    
                    %the fixed-increment algorithm has underestimated I_shar
                    check=0;
                    % thus, try another version of the algorithm where the increment is optimized with a line search
                    line_search=1;
                    % reset the algorithm to the initial starting point
                    q=p;
                    %reset the iteration count
                    iter=1;
                
                end
                
            else
                
                if line_search==0
                    %fixed increment
                    gamma_k=2/(iter+2);% this is the simplest version of Franke-Wolf algorithm.
                    
                else
                    % increment optimized with a line-search
                    
                    gamma=0;
                    gamma_k=0;
                    q(q<0)=0;
                    
                    q12 = sum(q, 3);
                    I_xy = q12 .* log2(q12 ./ repmat(sum(q12), [dimx 1]) ./ repmat(sum(q12,2), [1 dimy]));
                    I_xy = sum(I_xy(q12 > 0));
                    
                    %       I_q(X:Y|Z)
                    I_cond_xy_z = q .* log2(q ./ repmat(sum(sum(q), 2), [dimx dimy 1]) ./ ...
                        ( repmat(sum(q,2), [1 dimy 1]) ./ repmat(sum(sum(q), 2), [dimx dimy 1]) .* ...
                        repmat(sum(q),   [dimx 1 1]) ./ repmat(sum(sum(q), 2), [dimx dimy 1]) ) );
                    I_cond_xy_z = sum(I_cond_xy_z(q > 0));
                    
                    co_I_prev=I_xy-I_cond_xy_z;
                    
                    while gamma<1 %gamma=1 just gives you p_k again
                        
                        q_search=q+gamma*(p_k-q);
                        q_search(q_search<0)=0;
                        
                        q12 = sum(q_search, 3);
                        I_xy = q12 .* log2(q12 ./ repmat(sum(q12), [dimx 1]) ./ repmat(sum(q12,2), [1 dimy]));
                        I_xy = sum(I_xy(q12 > 0));
                        
                        %       I_q(X:Y|Z)
                        I_cond_xy_z = q_search .* log2(q_search ./ repmat(sum(sum(q_search), 2), [dimx dimy 1]) ./ ...
                            ( repmat(sum(q_search,2), [1 dimy 1]) ./ repmat(sum(sum(q_search), 2), [dimx dimy 1]) .* ...
                            repmat(sum(q_search),   [dimx 1 1]) ./ repmat(sum(sum(q_search), 2), [dimx dimy 1]) ) );
                        I_cond_xy_z = sum(I_cond_xy_z(q_search > 0));
                        
                        co_I_search=I_xy-I_cond_xy_z;
                        
                        if  co_I_search>co_I_prev
                            
                            gamma_k=gamma;
                            gamma=1;
                            co_I_prev=co_I_search;
                        end
                        
                        gamma=gamma+0.01;
                        
                    end
                    
                end
                
                q=q+gamma_k*(p_k-q); % update the q for next iteration
                
                coeff_prev=coeff_prev+gamma_k*(coeff_tot-coeff_prev); % update coeff_prev for next iteration
                
                iter=iter+1;
                
            end
            
        end
        
        I_shar=co_I; % get the output redundancy from the last coI_q
        
        
        I_xz=0;
        for i=1:dimx
            for j=1:dimz
                
                if sum(p(i,:,j))>0
                    I_xz=I_xz+sum(p(i,:,j))*log2(sum(p(i,:,j))/(sum(sum(p(:,:,j)))*sum(sum(p(i,:,:)))));
                end
                
            end
        end
        
        I_unx=I_xz-I_shar;% first output unique
        
        I_yz=0;
        for k=1:dimy
            for j=1:dimz
                
                if sum(p(:,k,j))>0
                    I_yz=I_yz+sum(p(:,k,j))*log2(sum(p(:,k,j))/(sum(sum(p(:,:,j)))*sum(sum(p(:,k,:)))));
                end
                
            end
        end
        
        
        I_uny=I_yz-I_shar;%second output unique
        
        I_xy_z=0;
        for i=1:dimx
            for k=1:dimy
                for j=1:dimz
                    
                    if p(i,k,j)>0
                        I_xy_z=I_xy_z+p(i,k,j)*log2(p(i,k,j)/(sum(p(i,k,:))*sum(sum(p(:,:,j)))));
                    end
                    
                end
            end
        end
        
        I_syn=I_xy_z-I_xz-I_yz+I_shar;%output synergy
        
    end
    
end
