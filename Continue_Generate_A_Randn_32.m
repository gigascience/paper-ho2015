clear A00

nmax = 8000;
p = 10000;

A00 = zeros(nmax,p);
A00 = randn(nmax,p);

random0 = zeros(nmax,1);
random0 = randn(nmax,1);

dlmwrite('A_Matrix_con',A00,'delimiter',',')
dlmwrite('RandomVector_con',random0,'delimiter',',')
