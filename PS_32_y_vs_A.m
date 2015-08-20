% This program is to plot y against A in the Promiscuous (PS) model with 
% s=3 and s'=2, using the synthetic data.

h2 = 0.7;
s = 3;
sprime = 2;
vare = 0.3;
p = 10000;
interval = 200;
num_min = 1;
num_max = 5;

load A_Matrix
load RandomVector
nmax = 8000;

nshuffle_max = 100;

ResidualVar = zeros(nshuffle_max,1);
nstar = zeros(nshuffle_max,1);
GXvar = zeros(nshuffle_max,1);
GXL2 = zeros(nshuffle_max,1);
xhat_first_collect0 = zeros(nshuffle_max*(s+sprime)+nshuffle_max,1);
xhat_first_collect = zeros(nshuffle_max*(s+sprime)+nshuffle_max,1);
error_length0 = 0;

for nshuffle = 1:nshuffle_max
        
 nshuffle   
      
 NLcoeff = zeros(sprime+sprime/2,1);
 NLcoeff_max = 1.5;
 NLcoeff_min = 0.5;
 NLcoeff = (NLcoeff_max - NLcoeff_min)*rand(sprime+sprime/2,1) + NLcoeff_min; 
 
 genotype = zeros(nmax,p);
 
 genotype = A_Matrix(:,randperm(p));
   
    n = nmax;

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
random(1:n) = RandomVector(1:n); 
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

if nnz(xhat_first0(4)) == 0 
    
    xhat_first0
    
    AAA = A(:,4); 

break

end

end

x1 = AAA;
y1 = y;

plot(x1,y1,'ro')
set(gca,'FontSize',11)
xlabel('standardized locus value','interpreter','latex','FontSize',20)
ylabel('y','interpreter','latex','FontSize',20)
h = lsline;
set(h,'color','b')
get(0,'screensize')
saveas(gcf,'PS_y_vs_locus.jpg')









