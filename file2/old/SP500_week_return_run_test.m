load('Alpha.mat')
load('SP500_return.mat');
load('SP500_return_cov.mat');
load('SP500_return_mean.mat');
load('variable_select.mat');

p0=250;n0=500;rep1=10;rep2=1000;sigma0=1;y=p0/n0;
V_sel=variable_select(1:rep1,1:p0);

i=7;
    variable_sel=V_sel(i,:);
    or_lrdata=SP500_return(:,variable_sel);
    mean=(SP500_return_mean(1,variable_sel))';
    S=SP500_return_cov(variable_sel,variable_sel);inv_S=inv(S);
    Cor=S./sqrt(diag(S)*diag(S)');
   % [Eigenvector_S,Eigenvalues_S]=eig(S);
    [Eigenvector_cor,Eigenvalues_cor]=eig(Cor);
    eigen_cor=diag(Eigenvalues_cor);
    Max=max(max(Eigenvalues_cor));
    
    par_1_invS_mu=ones(1,p0)*inv_S*mean;
    par_mu_invS_mu=mean'*inv_S*mean;
    par_1_invS_1=ones(1,p0)*inv_S*ones(p0,1);
    
    cond_plug=sigma0*par_1_invS_mu/sqrt(par_mu_invS_mu);
    b_plug=sqrt((par_1_invS_1*sigma0-1)/(par_mu_invS_mu*par_1_invS_1-(par_1_invS_mu)^2));
    
   num_plug=(cond_plug<1);
    
    if num_plug==1;
        x_plug=1; 
        return_plug=sqrt(par_mu_invS_mu);
    else x_plug=-1;
        return_plug=par_1_invS_mu/par_1_invS_1+b_plug*(par_mu_invS_mu-par_1_invS_mu^2/par_1_invS_1);   
    end
   %plug_return is end
   
   %spectral-corrected 
   alpha=Alpha(1,i);global alpha
   if Max>Eigenvalues_cor(1,1);
      F_vec=((1:p0)-0.1)/p0; x=1;
   else
   F_vec=1-((1:p0)-0.1)/p0;x=0;
   end
   est_cr=qutile(F_vec);
   est_cr_M=diag(est_cr');
   
   Est_cr=Eigenvector_cor*est_cr_M*Eigenvector_cor';
   S_spectral=sqrt(diag(S));
   S_sp=Est_cr.*(S_spectral*S_spectral');
   
   par_1_invSsp_1=ones(1,p0)*inv(S_sp)*ones(p0,1);
   par_1_invSsp_mu=ones(1,p0)*inv(S_sp)*mean;
   par_mu_invSsp_mu=mean'*inv(S_sp)*mean;
   
   cond_sp=sigma0*par_1_invSsp_mu/sqrt(par_mu_invSsp_mu);
   b_sp=sqrt((par_1_invSsp_1*sigma0-1)/(par_mu_invSsp_mu*par_1_invSsp_1-(par_1_invSsp_mu)^2));
   
   num_sp=(cond_sp<1);
    
    if num_sp==1;
        x_sp=1; 
        return_sp=sqrt(par_mu_invSsp_mu);
    else x_sp=-1;
        return_sp=par_1_invSsp_mu/par_1_invSsp_1+b_sp*(par_mu_invSsp_mu-par_1_invSsp_mu^2/par_1_invSsp_1);   
    end
   %end for spectral corrected method