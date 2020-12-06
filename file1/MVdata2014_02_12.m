
%第一步先输入数据DataMatrix_del_na，在dataMatrix中行代表不同的股票，列代表不同时间点。

LRn=size(DataMatrix);
ln=LRn(1,1);Rn=LRn(1,2);

%现在整理dataMatrix为logreturn形式的数据，lrdata means log return
lrdataup=DataMatrix(2:(ln-1),:);
lrdatadown=DataMatrix(3:(ln),:);
lrdata=log(lrdatadown)-log(lrdataup);
NP=(size(lrdata));N=NP(1,1);p=NP(1,2); y=p/N;global y;

%下面计算样本些方差矩阵for lrdata
sample_mean=0;
for i=1:p;
    sample_mean(i)=sum(lrdata(:,i));
end
sample_mean=sample_mean'/N;
lrdata=lrdata';
Mean_matrix=sample_mean*((1:1:N)*0+1);

%S=(lrdata-Mean_matrix)*(lrdata-Mean_matrix)'/N;
S=(lrdata)*(lrdata)'/N;
diag_S=sqrt(diag(S)*(diag(S))');

correlation_S=S./diag_S;

[Eigenvector,Eigenvalues]=eig(S);
[Eigen_cr_vector,Eigen_cr_values]=eig(correlation_S);
eigenY=diag(Eigenvalues);global eigenY;

%eigenY1=eigenY*100;
%eigenY1=eigen(1:158);
%eigenX=0:10:100;

plot(1:p,eigenY,'.')
x=0:0.00001:0.2;
hist(eigenY(3:p,1),x);figure(gcf);







