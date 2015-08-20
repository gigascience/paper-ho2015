load Con_32_xhatfirst0
load Con_32_xhatfirst
load Con_32_0var
load Con_32_nstar
load Con_32_var
load Con_32_L2

first01 = Con_32_xhatfirst0;
first1 = Con_32_xhatfirst;
var01 = Con_32_0var - 0.3*ones(100,1);
nstar1 = Con_32_nstar;
var1 = Con_32_var - 0.3*ones(100,1);
Lerr1 = Con_32_L2;

smax1 = 5;

A1 = zeros(smax1,100);
B1 = zeros(smax1,100);

for j = 1:100
    
    A1(:,j) = first01((j-1)*smax1+j+1:(j-1)*smax1+j+smax1);     
    B1(:,j) = first1((j-1)*smax1+j+1:(j-1)*smax1+j+smax1);     
    
end    

x01 = zeros(100,1);
x1 = zeros(100,1);

for j = 1:100

    x01(j) = (smax1 - nnz(A1(:,j)))/smax1;
    x1(j) = (smax1 - nnz(B1(:,j)))/smax1;
    
end

plot(var01,var1,'ro')
set(gca,'FontSize',11)
xlabel('nonlinear variance $\sigma_{\rm NL}^2$','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.5 0 0.2])
get(0,'screensize')
saveas(gcf,'Con_var_vs_var0.jpg')

plot(var01,nstar1/5,'ro')
set(gca,'FontSize',11)
xlabel('nonlinear variance $\sigma_{\rm NL}^2$','interpreter','latex','FontSize',20)
ylabel('$n_\ast / (s+s^{\prime})$','interpreter','latex','FontSize',20)
axis([0 0.5 100 250])
get(0,'screensize')
saveas(gcf,'Con_nstar_vs_var0.jpg')

plot(x1,var1,'ro')
set(gca,'FontSize',11)
xlabel('fraction of step 1 zeros','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.3 0 0.2])
get(0,'screensize')
saveas(gcf,'Con_var_vs_zero.jpg')

plot(x01,var1,'ro')
set(gca,'FontSize',11)
xlabel('fraction of model zeros','interpreter','latex','FontSize',20)
ylabel('residual variance  $\sigma_{\rm R}^2$','interpreter','latex','FontSize',20)
axis([0 0.05 0 0.2])
get(0,'screensize')
saveas(gcf,'Con_var_vs_0zero.jpg')
