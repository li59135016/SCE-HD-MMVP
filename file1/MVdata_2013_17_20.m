%��һ������������DataMatrix����dataMatrix���д���ͬ�Ĺ�Ʊ���д���ͬʱ��㡣

LRn=size(DataMatrix);
ln=LRn(1,1);Rn=LRn(1,2);

%��������dataMatrixΪlogreturn��ʽ�����ݣ�lrdata means log return
lrdataup=DataMatrix(2:(ln-1),:);
lrdatadown=DataMatrix(3:(ln),:);
lrdata=log(lrdatadown)-log(lrdataup);
NP=(size(lrdata));N=NP(1,1);p=NP(1,2);

%�����������Щ�������for lrdata
sample_mean=0;
for i=1:p;
    sample_mean(i)=sum(lrdata(:,i))/N;
end
sample_mean=sample_mean';
lrdata=lrdata';
Mean_matrix=sample_mean*((1:1:N)*0+1);

S=(lrdata-Mean_matrix)*(lrdata-Mean_matrix)'/N;

[Eigenvector,Eigenvalues]=eig(S);
eigenY=diag(Eigenvalues);
%eigenY1=eigenY*100;
%eigenY1=eigen(1:158);
%eigenX=0:10:100;

plot(1:p,eigenY,'.')
x=0:0.001:0.2;
hist(eigenY(3:p,1),x);%figure(gcf);
hist(lrdata(1,:),x);




