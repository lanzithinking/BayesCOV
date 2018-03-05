% This is to plot marginal posterior densities of covariances
% Prior is induced from squared-Dirichlet distribution

clear;
addpath('~/Documents/MATLAB/tight_subplot/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% parameters for setting
alpha=0.5;
sqdim=3; dim=sqdim*(sqdim+1)/2;

% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 600]);
nrow=floor(sqrt(dim)); ncol=ceil(dim/nrow);
ha=tight_subplot(nrow,ncol,[0.1,.05],[.09,.06],[.05,.04]);
lty={'-','--','-.',':'};


% plot estimated posteriors using Dlt-Spherical HMC
% load data
files = dir('./summary_sqdir');
nfiles = length(files) - 2;

found=false;
for j=1:nfiles
    if ~isempty(strfind(files(j+2).name,['_dim',num2str(dim),'_alpha',num2str(alpha),'_']))
        load(strcat('./summary_sqdir/', files(j+2).name));
        fprintf('%s loaded.\n',files(j+2).name);
        found=true;
    end
end
if found
    % convert posterior samples and calculate estimates
    N_postsamps=Nsamp;
    vecSigma_postsamp=zeros(N_postsamps,dim);
    sigma_sum=zeros(sqdim,1); Rho_sum=zeros(sqdim); Sigma_sum=zeros(sqdim);
    for s=1:N_postsamps
        [sigma_s,vecRho_s]=cnvrt_para(samp(s,1:sqdim)',samp(s,sqdim+1:end));
        sigma_sum=sigma_sum+sigma_s;
        Rho_s=ivech(vecRho_s,'row'); Rho_s=Rho_s+tril(Rho_s,-1)';
        Rho_sum=Rho_sum+Rho_s;
        Sigma_s=sigma_s.*Rho_s.*sigma_s';
        Sigma_sum=Sigma_sum+Sigma_s;
        vecSigma_postsamp(s,:)=vech(Sigma_s,'row');
    end
    sigma_est=sigma_sum./N_postsamps; Rho_est=Rho_sum./N_postsamps; Sigma_est=Sigma_sum./N_postsamps;
    % obtain sample correlation coefficients
    Sigma_mle=cov(data.y);
    for j=1:dim
        subplot(ha(j));
        [f,xi]=ksdensity(vecSigma_postsamp(:,j));
        plot(xi,f,'linestyle',lty{2},'linewidth',1.5,'displayname','Delta-SphHMC'); hold on;
        set(gca,'fontsize',14);
        [I,J]=ind_vech2sub(sqdim,j,'row');
        xlabel(['\sigma_{',num2str([I,J]),'}'],'fontsize',16);ylabel('Density','fontsize',16);
        plot(data.Sigma0(I,J).*ones(1,2),ylim,lty{1},'linewidth',2,'displayname','Truth'); hold on;
        plot(Sigma_mle(I,J).*ones(1,2),ylim,lty{3},'linewidth',1.5,'displayname','MLE'); hold on;
        h_lgd=legend('show','location','best');
        set(h_lgd,'fontsize',14);
    end
end
% save estimates
save('./summary_sqdir/post_recovertruth.mat','sigma_est','Rho_est','Sigma_est');

fig.PaperPositionMode = 'auto';
print(fig,'./summary_sqdir/post_recovertruth','-dpng','-r0');