load BD_5_false
load BD_50_false
load BD_100_false

load PS_32_false
load PS_3020_false
load PS_6040_false

load Data_5_false
load Data_32_false
load Con_32_false

A1 = BD_5_false;
A2 = BD_50_false;
A3 = BD_100_false;

A4 = PS_32_false;
A5 = PS_3020_false;
A6 = PS_6040_false;

A7 = Data_5_false;
A8 = Data_32_false;
A9 = Con_32_false;

F1 = (10000-5)*A1;
F1 = F1(F1>=0);
F1 = round(F1);
L1 = length(F1);

F2 = (25000-50)*A2;
F2 = F2(F2>=0);
F2 = round(F2);
L2 = length(F2);

F3 = (40000-100)*A3;
F3 = F3(F3>=0);
F3 = round(F3);
L3 = length(F3);

F4 = (10000-5)*A4;
F4 = F4(F4>=0);
F4 = round(F4);
L4 = length(F4);

F5 = (20000-50)*A5;
F5 = F5(F5>=0);
F5 = round(F5);
L5 = length(F5);

F6 = (30000-100)*A6;
F6 = F6(F6>=0);
F6 = round(F6);
L6 = length(F6);

F7 = (10000-5)*A7;
F7 = F7(F7>=0);
F7 = round(F7);
L7 = length(F7);

F8 = (10000-5)*A8;
F8 = F8(F8>=0);
F8 = round(F8);
L8 = length(F8);

F9 = (10000-5)*A9;
F9 = F9(F9>=0);
F9 = round(F9);
L9 = length(F9);

column_max = 8;

table = zeros(9,column_max);

for i = 1:9 
   
    for j = 1:column_max

if i==1    
    
table(i,j)=length(F1(F1==j-1))/L1;      

elseif i==2    
    
table(i,j)=length(F2(F2==j-1))/L2;      

elseif i==3    
    
table(i,j)=length(F3(F3==j-1))/L3;      

elseif i==4    
    
table(i,j)=length(F4(F4==j-1))/L4;      

elseif i==5    
    
table(i,j)=length(F5(F5==j-1))/L5;      

elseif i==6    
    
table(i,j)=length(F6(F6==j-1))/L6;      

elseif i==7    
    
table(i,j)=length(F7(F7==j-1))/L7;      

elseif i==8    
    
table(i,j)=length(F8(F8==j-1))/L8;      

else
    
table(i,j)=length(F9(F9==j-1))/L9;      

end

    end
    
end
  

dlmwrite('FPR_table',table,'delimiter',',')











