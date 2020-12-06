%第一步先输入数据DataMatrix，在dataMatrix中行代表不同的股票，列代表不同时间点。

LRn=size(DataMatrix);
ln=LRn(1,1);Rn=LRn(1,2);

sigma_0=1;%最初的风险要求
rep=10000;%boottrap重复的次数

%现在整理dataMatrix为logreturn形式的数据，lrdata means log return
lrdataup=DataMatrix(2:(ln-1),:);
lrdatadown=DataMatrix(3:(ln),:);
Lrdata=log(lrdatadown)-log(lrdataup);
lrdata=Lrdata(1:100,1:50);or_lrdata=lrdata;
NP=(size(lrdata));N=NP(1,1);p=NP(1,2);y=p/N;

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

%假设总体相关系数矩阵谱极限为 inverse cubic density. alpha=0.515805
alpha=0.763227;global alpha;
F_vec=1-((1:p)-0.1)/p;
est_cr=qutile(F_vec);%总体相关矩阵谱分布
est_cr_M=diag(est_cr');
Est_cr=Eigenvector_cor*est_cr_M*Eigenvector_cor';
S_spectral=sqrt(diag(S));
Est_cov=Est_cr.*(S_spectral*S_spectral');
[Eigenvector_Est_cov,Eigenvalues_Est_cov]=eig(Est_cov);

%现在计算plug-in return and allocation
par_1_S_mu=ones(1,p)*inv(S)*sample_mean;
par_1_S_1=ones(1,p)*inv(S)*ones(p,1);
par_mu_S_mu=sample_mean'*inv(S)*sample_mean;

b_S=sqrt((par_1_S_1-1)/(par_1_S_1*par_mu_S_mu-(par_1_S_mu)^2));
R_plug=par_1_S_mu/par_1_S_1+b_S*(par_mu_S_mu-(par_1_S_mu)^2/par_1_S_1);
c_S=inv(S)*ones(p,1)/par_1_S_1+b_S*(inv(S)*sample_mean-par_1_S_mu/par_1_S_1*inv(S)*ones(p,1));
risk_plug_Spectral=c_S'*inv(Est_cov)*c_S;
risk_plug_S=c_S'*inv(S)*c_S;

%现在计算spectral corrected estimation for return and risk;
par_1_Est_mu=ones(1,p)*inv(Est_cov)*sample_mean;
par_1_Est_1=ones(1,p)*inv(Est_cov)*ones(p,1);
par_mu_Est_mu=sample_mean'*inv(Est_cov)*sample_mean;

b_Est=sqrt((par_1_Est_1-1)/(par_1_Est_1*par_mu_Est_mu-(par_1_Est_mu)^2));
R_Est=par_1_Est_mu/par_1_Est_1+b_Est*(par_mu_Est_mu-(par_1_Est_mu)^2/par_1_Est_1);
c_Est=inv(Est_cov)*ones(p,1)/par_1_Est_1+b_Est*(inv(Est_cov)*sample_mean-par_1_Est_mu/par_1_Est_1*inv(Est_cov)*ones(p,1));
risk_Est_Spectral=c_Est'*inv(Est_cov)*c_Est;
risk_Est_S=c_Est'*inv(S)*c_Est;

%现在计算bootstrap estimation for return and risk;

R_plug_star=0;c_S_star=0;
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

R_boot=R_plug+(R_plug-R_plug_star)/sqrt(y);
c_boot=c_S+(c_S-c_S_star)/sqrt(y);
risk_boot_Spectral=c_boot'*Est_cov*c_boot;
risk_boot_S=c_boot'*Sample_cov*c_boot;

















