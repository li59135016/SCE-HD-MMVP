x=[1:13];
ans=combntns(x,5);
a=ans*34;
b=a-19;

n=length(a);
variable_select=[1:100];
for i=1:n
   v_s=[b(i,1):a(i,1),b(i,2):a(i,2),b(i,3):a(i,3),b(i,4):a(i,4),b(i,5):a(i,5)];
   variable_select=[variable_select;v_s] ;
   i
end

csvwrite('/Users/li59135016/Documents/MATLAB/variable_select.csv',variable_select);