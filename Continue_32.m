% This program is to study the Promiscuous (PS) model with s=3 and s'=2, 
% using the continuous synthetic data.

h2 = 0.7;
s = 3;
sprime = 2;
vare = 0.3;
p = 10000;
interval = 200;
num_min = 1;
num_max = 5;

load A_Matrix_con
load RandomVector_con
nmax = 8000;

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
         
 NLcoeff = zeros(sprime+sprime/2,1);
 NLcoeff_max = 1.5;
 NLcoeff_min = 0.5;
 NLcoeff = (NLcoeff_max - NLcoeff_min)*rand(sprime+sprime/2,1) + NLcoeff_min; 
 
 genotype = zeros(nmax,p);
 
 genotype = A_Matrix_con(:,randperm(p));
    
    n = nmax;

    A0 = zeros(n,p);

for i = 1:n

  A0(i,:) = abs(genotype(i,:));

end

A = zeros(n,p);
A = zscore(A0);

x = zeros(p,1);
x(1:s) = 1;
x(1:2:end) = -x(1:2:end);

y1 = A*x;
y2 = zeros(n,1);
for j = 1:n
    for k = 1:sprime
        if k < sprime/2+1
           y2(j) = y2(j)+ NLcoeff(k)*A(j,s+k)*A(j,s+k)+NLcoeff(sprime+k)*A(j,k)*A(j,s+k);
        end
        if k > sprime/2
           y2(j) = y2(j)+ NLcoeff(k)*A(j,s+k)*A(j,s+k); 
        end     
    end
end

random = zeros(n,1);
random(1:n) = RandomVector_con(1:n); 
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
     
 A0(i,:) = abs(genotype(i,:));

end
 
A = zeros(n,p);
A = zscore(A0);

x = zeros(p,1);
x(1:s) = 1;
x(1:2:end) = -x(1:2:end);

y1 = A*x;
y2 = zeros(n,1);
for j = 1:n
    for k = 1:sprime
        if k < sprime/2+1
           y2(j) = y2(j)+ NLcoeff(k)*A(j,s+k)*A(j,s+k)+NLcoeff(sprime+k)*A(j,k)*A(j,s+k);
        end
        if k > sprime/2
           y2(j) = y2(j)+ NLcoeff(k)*A(j,s+k)*A(j,s+k); 
        end     
    end
end

random = zeros(n,1);
random(1:n) = RandomVector_con(1:n); 
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

if Pvalue_n(num) <= 0.01*Pvalue_n(1) && abs(differP(num-1)) <= 0.01*abs(differP(1))
 
 ResidualVar(nshuffle) = sigmastar2; 
 nstar(nshuffle) = n;
 
 nshuffle
 n   
 ind;
  
 xhat_first = zeros(s+sprime,1);
 xhat_first(1:s+sprime) = xhat(1:s+sprime)
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

dlmwrite('Con_32_xhatfirst0',xhat_first_collect0,'delimiter',',')
dlmwrite('Con_32_xhatfirst',xhat_first_collect,'delimiter',',')
dlmwrite('Con_32_0var',ResidualVar,'delimiter',',')
dlmwrite('Con_32_nstar',nstar,'delimiter',',')
dlmwrite('Con_32_var',GXvar,'delimiter',',')
dlmwrite('Con_32_L2',GXL2,'delimiter',',')
dlmwrite('Con_32_false',false,'delimiter',',')

