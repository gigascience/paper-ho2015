load PS_32_xhatfirst0
load Data_32_xhatfirst0

load PS_32_xhatfirst
load Data_32_xhatfirst

load PS_32_0var
load Data_32_0var

load PS_32_nstar
load Data_32_nstar

load PS_32_var
load Data_32_var

load PS_32_L2
load Data_32_L2

first01 = PS_32_xhatfirst0;
first02 = Data_32_xhatfirst0;

first1 = PS_32_xhatfirst;
first2 = Data_32_xhatfirst;

% The penalization parameter used for real genome is "scaling_factor" times
% larger than that used for synthetic data:
scaling_factor = 2.0;

var01 = PS_32_0var - 0.3*ones(100,1);
var02 = Data_32_0var/scaling_factor - 0.3*ones(300,1);

nstar1 = PS_32_nstar;
nstar2 = Data_32_nstar;

nnz(nstar1);
nnz(nstar2);

var1 = PS_32_var - 0.3*ones(100,1);
var2 = Data_32_var - 0.3*ones(300,1);

Lerr1 = PS_32_L2;
Lerr2 = Data_32_L2;

smax1 = 5;
smax2 = 5;

A1 = zeros(smax1,100);
A2 = zeros(smax2,300);

B1 = zeros(smax1,100);
B2 = zeros(smax2,300);

for j = 1:100
    
    A1(:,j) = first01((j-1)*smax1+j+1:(j-1)*smax1+j+smax1);     
    B1(:,j) = first1((j-1)*smax1+j+1:(j-1)*smax1+j+smax1);     
    
end    

for j = 1:300
    
    A2(:,j) = first02((j-1)*smax2+j+1:(j-1)*smax2+j+smax2);     
    B2(:,j) = first2((j-1)*smax2+j+1:(j-1)*smax2+j+smax2);     
    
end    

x01 = zeros(100,1);
x02 = zeros(300,1);

x1 = zeros(100,1);
x2 = zeros(300,1);

for j = 1:100

    x01(j) = (smax1 - nnz(A1(:,j)))/smax1;
    x1(j) = (smax1 - nnz(B1(:,j)))/smax1;
    
end

for j = 1:300

    x02(j) = (smax2 - nnz(A2(:,j)))/smax2;
    x2(j) = (smax2 - nnz(B2(:,j)))/smax2;

end

plot(var01,var1,'ro',var02,var2,'b*')
set(gca,'FontSize',11)
xlabel('nonlinear variance $\sigma_{\rm NL}^2$','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.5 0 0.5])
get(0,'screensize')
saveas(gcf,'Data_PS_var_vs_var0.jpg')

plot(var01,nstar1/5,'ro',var02,nstar2/5,'b*')
set(gca,'FontSize',11)
xlabel('nonlinear variance $\sigma_{\rm NL}^2$','interpreter','latex','FontSize',20)
ylabel('$n_\ast / (s+s^{\prime})$','interpreter','latex','FontSize',20)
axis([0 0.5 100 250])
get(0,'screensize')
saveas(gcf,'Data_PS_nstar_vs_var0.jpg')

plot(x1,var1,'ro',x2,var2,'b*')
set(gca,'FontSize',11)
xlabel('fraction of step 1 zeros','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.5 0 0.5])
get(0,'screensize')
saveas(gcf,'Data_PS_var_vs_zero.jpg')

plot(x01,var1,'ro',x02,var2,'b*')
set(gca,'FontSize',11)
xlabel('fraction of model zeros','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.25 0 0.45])
get(0,'screensize')
saveas(gcf,'Data_PS_var_vs_0zero.jpg')
