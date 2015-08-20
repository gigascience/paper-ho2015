load BD_5_xhatfirst0
load BD_50_xhatfirst0
load BD_100_xhatfirst0

load BD_5_xhatfirst
load BD_50_xhatfirst
load BD_100_xhatfirst

load BD_5_0var
load BD_50_0var
load BD_100_0var

load BD_5_nstar
load BD_50_nstar
load BD_100_nstar

load BD_5_var
load BD_50_var
load BD_100_var

load BD_5_L2
load BD_50_L2
load BD_100_L2

first01 = BD_5_xhatfirst0;
first02 = BD_50_xhatfirst0;
first03 = BD_100_xhatfirst0;

first1 = BD_5_xhatfirst;
first2 = BD_50_xhatfirst;
first3 = BD_100_xhatfirst;

var01 = BD_5_0var - 0.3*ones(200,1);
var02 = BD_50_0var - 0.3*ones(100,1);
var03 = BD_100_0var - 0.3*ones(100,1);

nstar1 = BD_5_nstar;
nstar2 = BD_50_nstar;
nstar3 = BD_100_nstar;

nnz(nstar1);
nnz(nstar2);
nnz(nstar3);

var1 = BD_5_var - 0.3*ones(200,1);
var2 = BD_50_var - 0.3*ones(100,1);
var3 = BD_100_var - 0.3*ones(100,1);

Lerr1 = BD_5_L2;
Lerr2 = BD_50_L2;
Lerr3 = BD_100_L2;

smax1 = 5;
smax2 = 50;
smax3 = 100;

A1 = zeros(smax1,200);
A2 = zeros(smax2,100);
A3 = zeros(smax3,100);

B1 = zeros(smax1,200);
B2 = zeros(smax2,100);
B3 = zeros(smax3,100);

for j = 1:200
    
    A1(:,j) = first01((j-1)*smax1+j+1:(j-1)*smax1+j+smax1);     
    B1(:,j) = first1((j-1)*smax1+j+1:(j-1)*smax1+j+smax1);     
    
end    

for j = 1:100
    
    A2(:,j) = first02((j-1)*smax2+j+1:(j-1)*smax2+j+smax2);     
    A3(:,j) = first03((j-1)*smax3+j+1:(j-1)*smax3+j+smax3);     
    B2(:,j) = first2((j-1)*smax2+j+1:(j-1)*smax2+j+smax2);     
    B3(:,j) = first3((j-1)*smax3+j+1:(j-1)*smax3+j+smax3);     
    
end    

x01 = zeros(200,1);
x02 = zeros(100,1);
x03 = zeros(100,1);
x1 = zeros(200,1);
x2 = zeros(100,1);
x3 = zeros(100,1);

for j = 1:200

    x01(j) = (smax1 - nnz(A1(:,j)))/smax1;
    x1(j) = (smax1 - nnz(B1(:,j)))/smax1;
    
end

for j = 1:100

    x02(j) = (smax2 - nnz(A2(:,j)))/smax2;
    x03(j) = (smax3 - nnz(A3(:,j)))/smax3;
    x2(j) = (smax2 - nnz(B2(:,j)))/smax2;
    x3(j) = (smax3 - nnz(B3(:,j)))/smax3;

end

plot(var01,var1,'ro',var02,var2,'b*',var03,var3,'g+')
set(gca,'FontSize',11)
xlabel('nonlinear variance $\sigma_{\rm NL}^2$','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.5 0 0.25])
get(0,'screensize')
saveas(gcf,'BD_var_vs_var0.jpg')

plot(var01,nstar1/5,'ro',var02,nstar2/50,'b*',var03,nstar3/100,'g+')
set(gca,'FontSize',11)
xlabel('nonlinear variance $\sigma_{\rm NL}^2$','interpreter','latex','FontSize',20)
ylabel('$n_\ast / s$','interpreter','latex','FontSize',20)
axis([0 0.5 100 350])
get(0,'screensize')
saveas(gcf,'BD_nstar_vs_var0.jpg')

plot(x1,var1,'ro',x2,var2,'b*',x3,var3,'g+')
set(gca,'FontSize',11)
xlabel('fraction of step 1 zeros','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.5 0 0.25])
get(0,'screensize')
saveas(gcf,'BD_var_vs_zero.jpg')

plot(x01,var1,'ro',x02,var2,'b*',x03,var3,'g+')
set(gca,'FontSize',11)
xlabel('fraction of model zeros','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.25 0 0.25])
get(0,'screensize')
saveas(gcf,'BD_var_vs_0zero.jpg')
