load PS_32_xhatfirst0
load PS_3020_xhatfirst0
load PS_6040_xhatfirst0

load PS_32_xhatfirst
load PS_3020_xhatfirst
load PS_6040_xhatfirst

load PS_32_0var
load PS_3020_0var
load PS_6040_0var

load PS_32_nstar
load PS_3020_nstar
load PS_6040_nstar

load PS_32_var
load PS_3020_var
load PS_6040_var

load PS_32_L2
load PS_3020_L2
load PS_6040_L2

first01 = PS_32_xhatfirst0;
first02 = PS_3020_xhatfirst0;
first03 = PS_6040_xhatfirst0;

first1 = PS_32_xhatfirst;
first2 = PS_3020_xhatfirst;
first3 = PS_6040_xhatfirst;

var01 = PS_32_0var - 0.3*ones(100,1);
var02 = PS_3020_0var - 0.3*ones(100,1);
var03 = PS_6040_0var - 0.3*ones(100,1);

nstar1 = PS_32_nstar;
nstar2 = PS_3020_nstar;
nstar3 = PS_6040_nstar;

nnz(nstar1);
nnz(nstar2);
nnz(nstar3);

var1 = PS_32_var - 0.3*ones(100,1);
var2 = PS_3020_var - 0.3*ones(100,1);
var3 = PS_6040_var - 0.3*ones(100,1);

Lerr1 = PS_32_L2;
Lerr2 = PS_3020_L2;
Lerr3 = PS_6040_L2;

smax1 = 5;
smax2 = 50;
smax3 = 100;

A1 = zeros(smax1,100);
A2 = zeros(smax2,100);
A3 = zeros(smax3,100);

B1 = zeros(smax1,100);
B2 = zeros(smax2,100);
B3 = zeros(smax3,100);

for j = 1:100
    
    A1(:,j) = first01((j-1)*smax1+j+1:(j-1)*smax1+j+smax1);     
    A2(:,j) = first02((j-1)*smax2+j+1:(j-1)*smax2+j+smax2);     
    A3(:,j) = first03((j-1)*smax3+j+1:(j-1)*smax3+j+smax3);     
    B1(:,j) = first1((j-1)*smax1+j+1:(j-1)*smax1+j+smax1);     
    B2(:,j) = first2((j-1)*smax2+j+1:(j-1)*smax2+j+smax2);     
    B3(:,j) = first3((j-1)*smax3+j+1:(j-1)*smax3+j+smax3);     
    
end    

x01 = zeros(100,1);
x02 = zeros(100,1);
x03 = zeros(100,1);
x1 = zeros(100,1);
x2 = zeros(100,1);
x3 = zeros(100,1);

for j = 1:100

    x01(j) = (smax1 - nnz(A1(:,j)))/smax1;
    x02(j) = (smax2 - nnz(A2(:,j)))/smax2;
    x03(j) = (smax3 - nnz(A3(:,j)))/smax3;
    x1(j) = (smax1 - nnz(B1(:,j)))/smax1;
    x2(j) = (smax2 - nnz(B2(:,j)))/smax2;
    x3(j) = (smax3 - nnz(B3(:,j)))/smax3;

end

plot(var01,var1,'ro',var02,var2,'b*',var03,var3,'g+')
set(gca,'FontSize',11)
xlabel('nonlinear variance $\sigma_{\rm NL}^2$','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.5 0 0.5])
get(0,'screensize')
saveas(gcf,'PS_var_vs_var0.jpg')

plot(var01,nstar1/5,'ro',var02,nstar2/50,'b*',var03,nstar3/100,'g+')
set(gca,'FontSize',11)
xlabel('nonlinear variance $\sigma_{\rm NL}^2$','interpreter','latex','FontSize',20)
ylabel('$n_\ast / (s+s^{\prime})$','interpreter','latex','FontSize',20)
axis([0 0.5 100 250])
get(0,'screensize')
saveas(gcf,'PS_nstar_vs_var0.jpg')

plot(x1,var1,'ro',x2,var2,'b*',x3,var3,'g+')
set(gca,'FontSize',11)
xlabel('fraction of step 1 zeros','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.5 0 0.5])
get(0,'screensize')
saveas(gcf,'PS_var_vs_zero.jpg')

plot(x01,var1,'ro',x02,var2,'b*',x03,var3,'g+')
set(gca,'FontSize',11)
xlabel('fraction of model zeros','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.25 0 0.45])
get(0,'screensize')
saveas(gcf,'PS_var_vs_0zero.jpg')
