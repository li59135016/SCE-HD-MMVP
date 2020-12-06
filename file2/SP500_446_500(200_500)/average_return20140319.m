%Ê×ÏÈÊäÈëSP500_return_cov.csv;SP500_return_mean.csv;variable_select.csv(v_sel);Alpha.csv¡£

load('Alpha.mat')
load('SP500_return.mat');
load('SP500_return_cov.mat');
load('SP500_return_mean.mat');
load('variable_select.mat');

p0=200;n0=500;rep1=1000;rep2=1000;sigma0=1;y=p0/n0;
V_sel=variable_select(1:rep1,1:p0);

for i=1:rep1
    
    variable_sel=V_sel(i,:);
    or_lrdata=SP500_return(:,variable_sel);
    mean=(SP500_return_mean(1,variable_sel))';
    average_return(i)=sum(mean)/length(mean);
  i
end

mean_average=sum(average_return')/1000;
stand_average=sqrt(sum((average_return').^2)/1000-(mean_average).^2);
