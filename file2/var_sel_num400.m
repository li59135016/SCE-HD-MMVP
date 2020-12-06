n1=14;n2=9;n3=6;%n1为446分为n1个小组，n2为从n1中任选n2个小组，n3为每个小组中选择的元素个数。
x=[1:n1];
ans=combntns(x,n2);
a=ans*n3;
b=(ans-1)*n3+1;

n=length(a);
variable_select=[];
for i=1:n
   v_s=[];
   for j=1:n2
   v_s0=[b(i,j):a(i,j)];
   v_s=[v_s,v_s0];    
   end    
   %v_s=[b(i,1):a(i,1),b(i,2):a(i,2),b(i,3):a(i,3),b(i,4):a(i,4),b(i,5):a(i,5),b(i,6):a(i,6),b(i,7):a(i,7)];
   variable_select=[variable_select;v_s] ;
   i
end
variable_select1=variable_select(1:1001,1:50);
variable_select2=ones(1001,1)*[51:400];
variable_select400=[variable_select1,variable_select2];
csvwrite('/Users/li59135016/Dropbox/Markowitz MV/MV-data analysis/Matlab/file2/SP500_446_500(400_500)/variable_select400.csv',variable_select400);

%nchoosek(n,m)  n中取m个组合数

