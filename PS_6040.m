% This program is to study the Promiscuous (PS) model with s=60 and s'=40, 
% using the synthetic data.

% Heritability:
h2 = 0.7;
% Sparsity:
s = 60;
% Sparsity for nonlinear interactions:
sprime = 40;

vare = 0.3;
p = 30000;
interval = 4000;
num_min = 1;
num_max = 5;

% Loading the synthetic genotype and noise:
load A_Matrix_bigger
load RandomVector_bigger
nmax = 28000;

nshuffle_max = 100;

ResidualVar = zeros(nshuffle_max,1);
nstar = zeros(nshuffle_max,1);
GXvar = zeros(nshuffle_max,1);
GXL2 = zeros(nshuffle_max,1);
xhat_first_collect0 = zeros(nshuffle_max*(s+sprime)+nshuffle_max,1);
xhat_first_collect = zeros(nshuffle_max*(s+sprime)+nshuffle_max,1);
error_length0 = 0;
% The factor of "-100" below is arbitrary. What we need at the end are the
% non-negative entries extracted from "false". A large negative number like
% "-100" sharply distinguishes from the relevant non-negative numbers if we 
% look at "false" directly.
false = -100*ones(nshuffle_max,1);


% This "nshuffle" loop iterates nshuffle_max = 100 different realizations 
% of the model:
for nshuffle = 1:nshuffle_max
 
 % Generating random coefficients for the nonlinear interaction terms:   
 NLcoeff = zeros(sprime+sprime/2,1);
 NLcoeff_max = 1.5;
 NLcoeff_min = 0.5;
 NLcoeff = (NLcoeff_max - NLcoeff_min)*rand(sprime+sprime/2,1) + NLcoeff_min; 
 
 genotype = zeros(nmax,p);
 
 % Shuffling the different columns of the genotype matrix:
 genotype = A_Matrix_bigger(:,randperm(p));
   
    n = nmax;

    A0 = zeros(n,p);


% Euqating each column of A0 to that of genotype:    
for i = 1:n
    
    A0(i,:) = genotype(i,:);

end

A = zeros(n,p);
% Standardizing the A matrix:
A = zscore(A0);

% Createing a vector x with all entries being 0:
x = zeros(p,1);
% Re-assigning the first s entries of x to be 1:
x(1:s) = 1;
% Making all the odd elements negative:
x(1:2:end) = -x(1:2:end);

% The linear term:
y1 = A*x;

% The nonlinear interactions:
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
% Defining the random noise with "RandomVector_bigger" loaded above:
random(1:n) = RandomVector_bigger(1:n); 
% Defining y3, taking into account of the heritability:
y3 = sqrt(1-h2)*random;

% Optional parameters:
a = 1;
b = 1;

% Redefining y1, y2, y3 such that they are properly normalized. 
% This normalization automatically makes var(y) ~ 1:
y1prime = a*y1*sqrt(h2)/sqrt(var(a*y1+b*y2));
y2prime = b*y2*sqrt(h2)/sqrt(var(a*y1+b*y2));
y3prime = y3/sqrt(var(random));
y = y1prime + y2prime + y3prime;

% Standardizing y (could be redundant but no hurt to make sure):
y_STD = zscore(y);

% Calling the program "lasso_sv" to obtain x^* (asymptotic linear 
% effects vector):
xhat = lasso_sv(A,y_STD,h2,vare);

xhat_first0 = zeros(s+sprime,1);
% Storing the first s+sprime entries:
xhat_first0(1:s+sprime) = xhat(1:s+sprime);
xhat_first_collect0((nshuffle-1)*(s+sprime)+nshuffle) = nshuffle;
% Storing the first s+sprime entries for all of the different runs:
xhat_first_collect0((nshuffle-1)*(s+sprime)+nshuffle+1:(nshuffle-1)*(s+sprime)+nshuffle+(s+sprime)) = xhat_first0;

% Defining the nonlinear variance:
sigmastar2 = var(y_STD-A*xhat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The whole point of the above codes is to determine the nonlinear
% variance. From now on, we will run the program to study the phase 
% transition, residual variance, n*, etc. Some steps will be the same 
% as above. We will only add explanatory comments to those new codes.

% Clearing the similar symbols that have already been used above:
clear A0 A num n x random y1 y2 y3 y1prime y2prime y3prime y
clear a b y_STD xhat 

step = zeros(num_max,1);
Pvalue_n = zeros(num_max,1);
differP = zeros(num_max-1,1);

% This "num" loop iterates "num_max" values of n with an "interval" defined 
% at the begining of the program:
for num = num_min:num_max
           
   n = num*interval;
    
   A0 = zeros(n,p);

for i = 1:n
    
   A0(i,:) = genotype(i,:);

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
random(1:n) = RandomVector_bigger(1:n); 
y3 = sqrt(1-h2)*random;

a = 1;
b = 1;

y1prime = a*y1*sqrt(h2)/sqrt(var(a*y1+b*y2));
y2prime = b*y2*sqrt(h2)/sqrt(var(a*y1+b*y2));
y3prime = y3/sqrt(var(random));
y = y1prime + y2prime + y3prime;

y_STD = zscore(y);

% Running lasso to obtain xhat with vare taken to be the nonlinear
% variance determined above:
xhat = lasso_sv(A,y_STD,h2,sigmastar2);

% Finding the non-zero entries in xhat and keeping track of their indices:
ind = find(xhat);

% Determining the total number of non-zero entries (total number of 
% genetic markers that have non-zero support):
ind_max = length(ind);

PVal = zeros(ind_max,1);

for ijk = 1:ind_max

% Calling the program "lse" to compute the p-values of all genetic markers 
% that have non-zero support:    
P = lse(y_STD,A(:,ind(ijk)));

% Storing the p-values:
PVal(ijk) = P;

end

step(num) = n;

% Taking the median of the p-values:
Pvalue_n(num) = median(PVal);

% In order to calculate the difference below, we need to start from 
% num >= 2:
if num >= 2 && num <= num_max

differP(num-1) = Pvalue_n(num) - Pvalue_n(num-1);

% Imposing that the median p-value and the absolute value of its first
% derivative are both 10^6 times smaller than the corresponding quantities
% when the scanning process first starts:
if Pvalue_n(num) <= 0.000001*Pvalue_n(1) && abs(differP(num-1)) <= 0.000001*abs(differP(1))
 
 ResidualVar(nshuffle) = sigmastar2; 
 nstar(nshuffle) = n;
 
 % Printing "nshuffle" and "n" to see the progress of the running:
 nshuffle
 n   
 ind;
 
 xhat_first = zeros(s+sprime,1);
 % Storing the first s+sprime entries:
 xhat_first(1:s+sprime) = xhat(1:s+sprime);
 xhat_first_collect((nshuffle-1)*(s+sprime)+nshuffle) = nshuffle;
 % Storing the first s+sprime entries for all of the different runs:
 xhat_first_collect((nshuffle-1)*(s+sprime)+nshuffle+1:(nshuffle-1)*(s+sprime)+nshuffle+(s+sprime)) = xhat_first;   

 % Calculating the false positive rates in each run. The distribution of
 % probabilities for the false positives will be determined in the program
 % "Table.m"
 ind_xhat_first = find(xhat_first);
 ind_max0 = length(ind_xhat_first);
 false(nshuffle) = (ind_max-ind_max0)/(p-(s+sprime));
 
 % Defining the total number of genetic markers that have non-zero support
 % to be "seff":
 seff = ind_max;
 G0 = zeros(n,2*seff+(seff^2-seff)/2); 

% The following "for" loop builds the effective genotype matrix formed 
% over the subspace spanned by those genetic markers that have non-zero 
% support:
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

% Standarding G0:
G = zscore(G0);

% Running lasso to obtain bigXhat0:
bigXhat0 = lasso_sv(G,y_STD,h2,vare);
% Setting vareGX = residual variance of G*X: 
vareGX = var(y_STD - G*bigXhat0);

% Running lasso again with vareGX = residual variance of G*X:
bigXhat = lasso_sv(G,y_STD,h2,vareGX);

% The value of BigY = Y should be the same as y. This has been
% double-cheked to ensure that the code is correct:
BigY = G*bigXhat + y3prime;

% Calulating the L2 errors:
L2error = norm(BigY-y_STD,2)/norm(y_STD,2);

% Storing the residual variance defined in the paper:
GXvar(nshuffle) = var(y_STD - G*bigXhat);
GXL2(nshuffle) = L2error;

% Breaking the loop caculations when the follwing two conditions are 
% satisfied: the median p-value and the absolute value of its first 
% derivative are both 10^6 times smaller than the corresponding quantities 
% when the scanning process first starts. This brings the calculation back 
% to the beginning of the "nshuffle" loop, and the next value of 
% "nshuffle" will be executed: 
break
    
end

end

end

end

% Writing the files:
dlmwrite('PS_6040_xhatfirst0',xhat_first_collect0,'delimiter',',')
dlmwrite('PS_6040_xhatfirst',xhat_first_collect,'delimiter',',')
dlmwrite('PS_6040_0var',ResidualVar,'delimiter',',')
dlmwrite('PS_6040_nstar',nstar,'delimiter',',')
dlmwrite('PS_6040_var',GXvar,'delimiter',',')
dlmwrite('PS_6040_L2',GXL2,'delimiter',',')
dlmwrite('PS_6040_false',false,'delimiter',',')

