%假设总体相关系数矩阵谱极限为 inverse cubic density. alpha=0.515805
alpha=0.276478; global alpha;
%F_vec=1-((1:p)-0.1)/p;
F_vec=((1:p)-0.1)/p;
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

y=p0/n0;

R_plug_star=R_plug_star/rep;
c_S_star=c_S_star/rep;

R_boot=R_plug+(R_plug-R_plug_star)/sqrt(y);
c_boot=c_S+(c_S-c_S_star)/sqrt(y);
risk_boot_Spectral=c_boot'*Est_cov*c_boot;
risk_boot_S=c_boot'*Sample_cov*c_boot;

Return=[R_plug,R_boot,R_Est]';
Risk_spectral=[risk_plug_Spectral,risk_boot_Spectral,risk_Est_Spectral]';

Return

Risk_spectral


