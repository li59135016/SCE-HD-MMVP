LRn=size(DataMatrix);
ln=LRn(1,1);Rn=LRn(1,2);

sigma_0=1;%

%现在整理dataMatrix为logreturn形式的数据，lrdata means log return
lrdataup=DataMatrix(2:(ln-1),:);
lrdatadown=DataMatrix(3:(ln),:);
or_lrdata=log(lrdatadown)-log(lrdataup);
NP=(size(or_lrdata));N=NP(1,1);p=NP(1,2);

rep=10;R_plug_star=0;c_S_star=0;
for j=1:rep
   loc=sort(fix(rand(1,N)*N)+1);
   lrdata=or_lrdata(loc,1:p);
   sample_mean=0;
for i=1:p;
    sample_mean(i)=sum(lrdata(:,i))/N;
end
sample_mean=sample_mean';
lrdata=lrdata';
Mean_matrix=sample_mean*((1:1:N)*0+1);

S=(lrdata-Mean_matrix)*(lrdata-Mean_matrix)'/N;
[Eigenvector,Eigenvalues]=eig(S);

par_1_S_mu=ones(1,p)*inv(S)*sample_mean;
par_1_S_1=ones(1,p)*inv(S)*ones(p,1);
par_mu_S_mu=sample_mean'*inv(S)*sample_mean;

b_S=sqrt((par_1_S_1-1)/(par_1_S_1*par_mu_S_mu-(par_1_S_mu)^2));
R_plug_star=R_plug_star+par_1_S_mu/par_1_S_1+b_S*(par_mu_S_mu-(par_1_S_mu)^2/par_1_S_1);
c_S_star=c_S_star+inv(S)*ones(p,1)/par_1_S_1+b_S*(inv(S)*sample_mean-par_1_S_mu/par_1_S_1*inv(S)*ones(p,1));
j
end
R_plug_star=R_plug_star/rep;
c_S_star=c_S_star/rep;

