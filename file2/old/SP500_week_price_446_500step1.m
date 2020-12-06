%先输入SP500_data_price_week_446_500.mat,
LRn=size(SP500weekprice466500);
ln=LRn(1,1);Rn=LRn(1,2);

%现在整理dataMatrix为logreturn形式的数据，lrdata means log return
lrdataup=SP500weekprice466500(1:(ln-1),:);
lrdatadown=SP500weekprice466500(2:(ln),:);
lrdata=log(lrdatadown)-log(lrdataup);
%计算log return

Varnumber=[1:Rn];
Varnumber=Varnumber(randperm(numel(Varnumber)));
Varname=SP500weekprice466500name(Varnumber);

p0=250;n0=500;rep=10000;
Lrdata=lrdata(:,Varnumber(1:p0));
Lrdataname=Varname(1:p0);
LrData=Lrdata;or_lrdata=Lrdata;%lrdata=Lrdata;or_lrdata=lrdata;

p=p0;N=n0;

sample_mean=0;
for i=1:p;
    sample_mean(i)=sum(LrData(:,i))/N;
end
sample_mean=sample_mean';
Lrdata=LrData';
Mean_matrix=sample_mean*((1:1:N)*0+1);

S=(Lrdata-Mean_matrix)*(Lrdata-Mean_matrix)'/N;
Sample_cov=S;
Cor=S./sqrt(diag(S)*diag(S)');

[Eigenvector,Eigenvalues]=eig(S);
[Eigenvector_cor,Eigenvalues_cor]=eig(Cor);
eigen_S=diag(Eigenvalues);
eigen_cor=diag(Eigenvalues_cor);
eigen_csv=eigen_cor';

 csvwrite('/Users/li59135016/Dropbox/Markowitz MV/MV-data analysis/Matlab/file2/eigen_csv.csv',eigen_csv);
 csvwrite('/Users/li59135016/Dropbox/Markowitz MV/MV-data analysis/Matlab/file2/SP500_return.csv',lrdata);

Max_S=max(eigen_S);Max_cor=max(eigen_cor);
x=0:0.01:Max_cor;
hist(eigen_cor,x)


