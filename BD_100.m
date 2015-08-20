% This program is to study the Block-Diagonal (BD) model with s=100, using
% the synthetic data.

clear A error TruncatedGeno TruncatedCorr
clear linearNum nonlinearNum linearVar nonlinearVar nstar

h2 = 0.7;
s = 100;
sprime = 0;
vare = 0.3;
p = 40000;
interval = 6000;
num_min = 1;
num_max = 5;

load A_Matrix_bigger_BD
load RandomVector_bigger_BD
nmax = 40000;

nshuffle_max = 100;

ResidualVar = zeros(nshuffle_max,1);
nstar = zeros(nshuffle_max,1);
GXvar = zeros(nshuffle_max,1);
GXL2 = zeros(nshuffle_max,1);
xhat_first_collect0 = zeros(nshuffle_max*(s+sprime)+nshuffle_max,1);
xhat_first_collect = zeros(nshuffle_max*(s+sprime)+nshuffle_max,1);
error_length0 = 0;
false = -100*ones(nshuffle_max,1);

for nshuffle = 1:nshuffle_max
        
 clear genotype
 clear A0
  
 genotype = zeros(nmax,p);
 
 genotype = A_Matrix_bigger_BD(:,randperm(p));
    
    n = nmax;

    A0 = zeros(n,p);

for i = 1:n
    
    A0(i,:) = genotype(i,:);

end

A = zeros(n,p);
A = zscore(A0);

aa = 1.5 + 0.5*randn(s,1);
bb = 1.0 + 0.2*randn(s,1);
cc = 0.5 + 0.1*randn(s,1);

aa(1:2:end) = -aa(1:2:end);
aa = aa(randperm(s));

bb(1:2:end) = -bb(1:2:end);
bb = bb(randperm(s));

cc(1:2:end) = -cc(1:2:end);
cc = cc(randperm(s));

y1 = zeros(n,1);

for j = 1:n

 for k = 1:s    
    
  y1(j) = y1(j) + aa(k)*A(j,k);    
 
 end
 
end

y2 =zeros(n,1);

for j = 1:n
    
    for k = 1:s
    
        if k < s
        
            y2(j) = y2(j)+ bb(k)*A(j,k)*A(j,k)+cc(k)*A(j,k)*A(j,k+1);
        
        end
        
        if k > s-1
        
            y2(j) = y2(j)+ bb(k)*A(j,k)*A(j,k); 
        
        end
        
    end
    
end

random = zeros(n,1);
random(1:n) = RandomVector_bigger_BD(1:n); 
y3 = sqrt(1-h2)*random;

a = 1;
b = 1;

y1prime = a*y1*sqrt(h2)/sqrt(var(a*y1+b*y2));
y2prime = b*y2*sqrt(h2)/sqrt(var(a*y1+b*y2));
y3prime = y3/sqrt(var(random));
y = y1prime + y2prime + y3prime;

y_STD = zscore(y);

xhat = lasso_sv(A,y_STD,h2,vare);

xhat_first0 = zeros(s+sprime,1);
xhat_first0(1:s+sprime) = xhat(1:s+sprime);
xhat_first_collect0((nshuffle-1)*(s+sprime)+nshuffle) = nshuffle;
xhat_first_collect0((nshuffle-1)*(s+sprime)+nshuffle+1:(nshuffle-1)*(s+sprime)+nshuffle+(s+sprime)) = xhat_first0;

sigmastar2 = var(y_STD-A*xhat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear A0 A num n x random y1 y2 y3 y1prime y2prime y3prime y
clear a b y_STD xhat 

step = zeros(num_max,1);
Pvalue_n = zeros(num_max,1);
differP = zeros(num_max-1,1);

for num = num_min:num_max
           
   n = num*interval;
    
   A0 = zeros(n,p);

for i = 1:n
    
   A0(i,:) = genotype(i,:);

end

A = zeros(n,p);
A = zscore(A0);

y1 = zeros(n,1);

for j = 1:n

 for k = 1:s    
    
  y1(j) = y1(j) + aa(k)*A(j,k);    
 
 end
 
end

y2 =zeros(n,1);

for j = 1:n
    
    for k = 1:s
    
        if k < s
        
            y2(j) = y2(j)+ bb(k)*A(j,k)*A(j,k)+cc(k)*A(j,k)*A(j,k+1);
        
        end
        
        if k > s-1
        
            y2(j) = y2(j)+ bb(k)*A(j,k)*A(j,k); 
        
        end
        
    end
    
end

random = zeros(n,1);
random(1:n) = RandomVector_bigger_BD(1:n); 
y3 = sqrt(1-h2)*random;

a = 1;
b = 1;

y1prime = a*y1*sqrt(h2)/sqrt(var(a*y1+b*y2));
y2prime = b*y2*sqrt(h2)/sqrt(var(a*y1+b*y2));
y3prime = y3/sqrt(var(random));
y = y1prime + y2prime + y3prime;

y_STD = zscore(y);
xhat = lasso_sv(A,y_STD,h2,sigmastar2);

ind = find(xhat);

ind_max = length(ind);

PVal = zeros(ind_max,1);

for ijk = 1:ind_max

P = lse(y_STD,A(:,ind(ijk)));

PVal(ijk) = P;

end

step(num) = n;

Pvalue_n(num) = median(PVal);

if num >= 2 && num <= num_max

differP(num-1) = Pvalue_n(num) - Pvalue_n(num-1);

if Pvalue_n(num) <= 0.000001*Pvalue_n(1) && abs(differP(num-1)) <= 0.000001*abs(differP(1))
 
 ResidualVar(nshuffle) = sigmastar2; 
 nstar(nshuffle) = n;
 
 nshuffle
 n   
 ind;
  
 xhat_first = zeros(s+sprime,1);
 xhat_first(1:s+sprime) = xhat(1:s+sprime);
 xhat_first_collect((nshuffle-1)*(s+sprime)+nshuffle) = nshuffle;
 xhat_first_collect((nshuffle-1)*(s+sprime)+nshuffle+1:(nshuffle-1)*(s+sprime)+nshuffle+(s+sprime)) = xhat_first;   
 
 ind_xhat_first = find(xhat_first);
 ind_max0 = length(ind_xhat_first);
 false(nshuffle) = (ind_max-ind_max0)/(p-(s+sprime));

 seff = ind_max;
 G0 = zeros(n,2*seff+(seff^2-seff)/2); 
  
for rr = 1:n

    for tt = 1:2*seff
    
      if tt < seff+1
        
        G0(rr,tt) = A(rr,ind(tt));
      
      end
      
      if tt > seff && tt < 2*seff+1
        
        G0(rr,tt) = A(rr,ind(tt-seff))*A(rr,ind(tt-seff));
        
      end    
     
    end
    
    column_index0 = 0;
    
    for uu = 1:seff-1    
        
       for vv = 1:seff-uu  
                   
         column_index = column_index0 + vv;
         
         G0(rr,2*seff+column_index) = A(rr,ind(uu))*A(rr,ind(uu+vv));
             
       end
       
         column_index0 = column_index;
              
     end
        
end

G = zscore(G0);

bigXhat0 = lasso_sv(G,y_STD,h2,vare);
vareGX = var(y_STD - G*bigXhat0);
bigXhat = lasso_sv(G,y_STD,h2,vareGX);

BigY = G*bigXhat + y3prime;

L2error = norm(BigY-y_STD,2)/norm(y_STD,2);

GXvar(nshuffle) = var(y_STD - G*bigXhat);
GXL2(nshuffle) = L2error;

break
    
end

end

end

end

dlmwrite('BD_100_xhatfirst0',xhat_first_collect0,'delimiter',',')
dlmwrite('BD_100_xhatfirst',xhat_first_collect,'delimiter',',')
dlmwrite('BD_100_0var',ResidualVar,'delimiter',',')
dlmwrite('BD_100_nstar',nstar,'delimiter',',')
dlmwrite('BD_100_var',GXvar,'delimiter',',')
dlmwrite('BD_100_L2',GXL2,'delimiter',',')
dlmwrite('BD_100_false',false,'delimiter',',')
