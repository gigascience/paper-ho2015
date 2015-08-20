clear A00

nmax = 18000;
p = 20000;

q_max = 0.49;
q_min = 0.05;
q = (q_max-q_min)*rand(p,1) + q_min;

A00 = zeros(nmax,p);

for i = 1:p
    % The number of entries with "0":
    n0 = round(nmax*(1-q(i))^2);
    % The number of entries with "1": 
    n1 = round(nmax*2*q(i)*(1-q(i)));
    % The number of entries with "2":
    n2 = nmax-n0-n1;
    % This assigns the entries "0", "1" and "2" into the i^th column of A:
    A00(1:n0,i) = 0;
    A00(n0+1:n0+n1,i) = 1;
    A00(n0+n1+1:nmax,i) = 2;
    % This is to shuffle the entries in the i^th column of A:
    A00(1:nmax,i) = A00(randperm(nmax),i);
end

random0 = zeros(nmax,1);
random0 = randn(nmax,1);

dlmwrite('A_Matrix_big',A00,'delimiter',',')
dlmwrite('RandomVector_big',random0,'delimiter',',')

