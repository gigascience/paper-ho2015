%References:
%1. Friedman,  J. H., Hastie, T., Hoefling H., and Tibshirani, R. Pathwise Coordinate Optimization. 
%Annals Applied Statistics 1, 302-332 (2007)
%2. Friedman,  J. H., Hastie, T. and Tibshirani, R. Regularized Paths for Generalized Linear Models 
%via Coordinate Descent. Journal of Statistical Software, 33(1) (2010)
%see section 2.5 of ref. 2 for notes on warm start (the lambda progression below) 

function xhat=lasso(A,y,h2,vare)
    [n,p]=size(A); %A is expected to be a matrix of standard normal columns, with subjects as 
                   %rows (n) and genetic variants as columns (p)
    A=A/sqrt(n);
    y=y/sqrt(n); %y is expected to be a standard normal column 
    %calc lambda max
    lambda_max=max(abs(A'*y));
    %calc lambda min
        %use the median maximum value of A transpose times a noise vector z;
        %median taken from 1,000 samples of z (1,000 was arbirtrarily chosen)
        %z is scaled by sige (and 1/sqrt(n) as above for A and y); also included is a
        %sampling error term 1/sqrt(n), found to be important when sige=0
    %vare=1-h2;
    %vare=0.0;
    z=sqrt(vare+1/n)*randn(n,1000)/sqrt(n);
    lambda_min=median(max(abs(A'*z)));    
    %sige=sqrt(1-h2);
    %z=(sige+1/sqrt(n))*randn(n,1000)/sqrt(n);
    %lambda_min=median(max(abs(A'*z)));
    %calc lambda series
    loghi = log(lambda_max);
    loglo = log(lambda_min);
    logrange = loghi - loglo;
    nLambda=100;
    interval = -logrange/(nLambda-1);
    lambda_ = exp(loghi:interval:loglo)';
    ITER=0; %global (master loop) iteration count, not used for other variables
    for lambi=1:length(lambda_)
        if lambi==1
            xhat=A'*y; %initialize coefficients by LSE if first lambda; otherwise, use fit from prior lambda
        end
        lambda=lambda_(lambi);
        r=y-A*xhat;
        kill=0; %kill variable for lasso subloop
        iter=0; %iteration count lasso subloop, used for indexing error
        err=[];
        %initialize active set, all
        active=ones(length(xhat),1); %active keeps track of current nonzeros; for iter==1, scan all parameters
        while kill==0
            iter=iter+1;
            ind=find(active~=0); %grab index of current nonzeros, and fit only these
            for ji=1:length(ind)
                j=ind(ji);
                xjold = xhat(j);
                xhat(j) = xjold+A(:,j)'*(r);
                xhat(j) = sign(xhat(j)).*max(abs(xhat(j))-lambda,0);
                r=r-A(:,j)*(xhat(j)-xjold); %update total residual by change in partial residual
                active(j)=xhat(j)~=0;
            end
            ind=find(xhat~=0);
            err(iter)=sum(r.^2)+lambda*sum(abs(xhat(ind)));
            if iter>1
               if 1-min([err(iter-1),err(iter)])/max([err(iter-1),err(iter)])<10^-4  | isfinite(err(iter))~=1 
                    %end lamba subloop if delta error is below threshold or if err==NaN
                  kill=1;
               end
            end
        end
        ITER=ITER+iter; %update global iteration count
    end
end