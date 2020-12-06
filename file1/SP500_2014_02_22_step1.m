sp500return446375%第一步先输入数据sp500_2014_02_22.mat，在dataMatrix中行代表不同的股票，列代表不同时间点。
p0=80;n0=400;
rep=100000;
Lrdata=SP50020140222(1:p0,1:n0);
Lrdata=Lrdata';
lrdata=Lrdata;or_lrdata=lrdata;
NP=(size(Lrdata));N=NP(1,1);p=NP(1,2);y=p/N;

%下面计算样本些方差矩阵for lrdata
sample_mean=0;
for i=1:p;
    sample_mean(i)=sum(lrdata(:,i))/N;
end
sample_mean=sample_mean';
lrdata=lrdata';
Mean_matrix=sample_mean*((1:1:N)*0+1);

S=(lrdata-Mean_matrix)*(lrdata-Mean_matrix)'/N;
Sample_cov=S;
Cor=S./sqrt(diag(S)*diag(S)');

[Eigenvector,Eigenvalues]=eig(S);
[Eigenvector_cor,Eigenvalues_cor]=eig(Cor);
eigen_S=diag(Eigenvalues);
eigen_cor=diag(Eigenvalues_cor);

%eigenY1=eigenY*100;
%eigenY1=eigen(1:158);
%eigenX=0:10:100;

%plot(1:p,eigen_cor,'.')
%x=0:0.05:38;
%hist(eigenY(3:p,1),x);
%hist(lrdata(1,:),x);
%hist(eigen_cor(1:p,1),x);figure(gcf)


















