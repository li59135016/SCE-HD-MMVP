%首先输入SP500_return_cov.csv;SP500_return_mean.csv;variable_select.csv(v_sel);Alpha.csv。

load('Alpha.mat')
load('SP500_return.mat');
load('SP500_return_cov.mat');
load('SP500_return_mean.mat');
load('variable_select.mat');

p0=100;n0=500;rep1=1000;rep2=1000;sigma0=1;y=p0/n0;
V_sel=variable_select(1:rep1,1:p0);



for i=1:rep1
    
    variable_sel=V_sel(i,:);
    or_lrdata=SP500_return(:,variable_sel);
    mean=(SP500_return_mean(1,variable_sel))';
    average_return(i)=sum(mean)/length(mean);
    S=SP500_return_cov(variable_sel,variable_sel);inv_S=inv(S);
    Cor=S./sqrt(diag(S)*diag(S)');
   % [Eigenvector_S,Eigenvalues_S]=eig(S);
    [Eigenvector_cor,Eigenvalues_cor]=eig(Cor);
    Max=max(max(Eigenvalues_cor));
    
    %plug_in estimation
    par_1_invS_mu(i)=ones(1,p0)*inv_S*mean;
    par_mu_invS_mu(i)=mean'*inv_S*mean;
    par_1_invS_1(i)=ones(1,p0)*inv_S*ones(p0,1);
    
    cond_plug(i)=sigma0*par_1_invS_mu(i)/sqrt(par_mu_invS_mu(i));
    b_plug(i)=sqrt((par_1_invS_1(i)*sigma0-1)/(par_mu_invS_mu(i)*par_1_invS_1(i)-(par_1_invS_mu(i))^2));
    
    num_plug(i)=(cond_plug(i)<1);
    
    if num_plug(i)==1;
        x_plug(i)=1; 
        return_plug(i)=sqrt(par_mu_invS_mu(i));
    else x_plug(i)=-1;
        return_plug(i)=par_1_invS_mu(i)/par_1_invS_1(i)+b_plug(i)*(par_mu_invS_mu(i)-par_1_invS_mu(i)^2/par_1_invS_1(i));   
    end
   %plug_return is end
   
   %spectral-corrected 
   alpha=Alpha(1,i);global alpha
   if Max>Eigenvalues_cor(1,1);
      F_vec=((1:p0)-0.1)/p0; x(i)=1;
   else
   F_vec=1-((1:p0)-0.1)/p0;x(i)=0;
   end
   est_cr=qutile(F_vec);
   est_cr_M=diag(est_cr');
   
   Est_cr=Eigenvector_cor*est_cr_M*Eigenvector_cor';
   S_spectral=sqrt(diag(S));
   S_sp=Est_cr.*(S_spectral*S_spectral');
   
   par_1_invSsp_1(i)=ones(1,p0)*inv(S_sp)*ones(p0,1);
   par_1_invSsp_mu(i)=ones(1,p0)*inv(S_sp)*mean;
   par_mu_invSsp_mu(i)=mean'*inv(S_sp)*mean;
   
   cond_sp(i)=sigma0*par_1_invSsp_mu(i)/sqrt(par_mu_invSsp_mu(i));
   b_sp(i)=sqrt((par_1_invSsp_1(i)*sigma0-1)/(par_mu_invSsp_mu(i)*par_1_invSsp_1(i)-(par_1_invSsp_mu(i))^2));
   
   num_sp(i)=(cond_sp(i)<1);
    
    if num_sp(i)==1;
        x_sp(i)=1; 
        return_sp(i)=sqrt(par_mu_invSsp_mu(i));
    else x_sp(i)=-1;
        return_sp(i)=par_1_invSsp_mu(i)/par_1_invSsp_1(i)+b_sp(i)*(par_mu_invSsp_mu(i)-par_1_invSsp_mu(i)^2/par_1_invSsp_1(i));   
    end
   %end for spectral corrected method
   
   %现在计算bootstrap estimation for return and risk;
   
   R_plug_star=0;c_S_star=0;N=n0;p=p0;
   
   for j=1:rep2
   loc=sort(fix(rand(1,N)*N)+1);
   lrdata=or_lrdata(loc,1:p);
   sample_mean=0;
   for t=1:p;
    sample_mean(t)=sum(lrdata(:,t))/N;
   end
   sample_mean=sample_mean';
   lrdata=lrdata';
   Mean_matrix=sample_mean*((1:1:N)*0+1);

   S=(lrdata-Mean_matrix)*(lrdata-Mean_matrix)'/N;
   [Eigenvector,Eigenvalues]=eig(S);

   par_1_S_mu=ones(1,p)*inv(S)*sample_mean;
   par_1_S_1=ones(1,p)*inv(S)*ones(p,1);
   par_mu_S_mu=sample_mean'*inv(S)*sample_mean;

   cond_star=sigma0*par_1_S_mu/sqrt(par_mu_S_mu);
   num_star=(cond_star<1);
   
   b_S=sqrt((par_1_S_1-1)/(par_1_S_1*par_mu_S_mu-(par_1_S_mu)^2));
   
   if num_star==1
       R_plug_star(j)=sqrt(par_mu_S_mu);
   else
   R_plug_star(j)=par_1_S_mu/par_1_S_1+b_S*(par_mu_S_mu-(par_1_S_mu)^2/par_1_S_1);
   end 
   
[i,j]
end
  return_plug_star(i)=sum(R_plug_star)/rep2;
  return_b(i)=return_plug(i)+(return_plug(i)-return_plug_star(i))/sqrt(y);
 
end

%return_combine=[return_plug;return_b;return_sp];
ans1=1:rep1;ans2=return_plug<return_sp;
return_combine=[ans1;ans2;return_plug;average_return;return_b;return_sp];
num_p_s=[num_plug;num_sp];
save('return_combine');
mean_combine=sum(return_combine')/1000;
stand_combine=sqrt(sum((return_combine').^2)/1000-(mean_combine).^2);
save('mean_combine');
save('stand_combine');