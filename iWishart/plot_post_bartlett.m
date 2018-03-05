% This is to plot marginal posterior densities of Gaussian-inverse-Wishart

clear;
addpath('~/Documents/MATLAB/tight_subplot/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% parameters for setting
sqdim=3; dim=sqdim*(sqdim+1)/2;
Psi=eye(sqdim); nu=sqdim;
N_prisamps=1e6;reweight=false;
if reweight
    ifrewt='_reweight';
else
    ifrewt='_noreweight';
end

% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 600]);
nrow=floor(sqrt(dim)); ncol=ceil(dim/nrow);
ha=tight_subplot(nrow,ncol,[0.1,.05],[.09,.06],[.05,.04]);
lty={'-','--','-.',':'};

% plot true posteriors vs estimated posteriors using Dlt-Spherical HMC
% load data
files = dir('./summary_bartlett');
nfiles = length(files) - 2;
found=false;
for i=1:nfiles
    if ~isempty(strfind(files(i+2).name,['_dim',num2str(dim),'_nu',num2str(nu),'_'])) && contains(files(i+2).name, ifrewt)
        load(strcat('./summary_bartlett/', files(i+2).name));
        fprintf('%s loaded.\n',files(i+2).name);
        found=true;
    end
end
if found
    % obtain data and sample from the true posterior
    Psi_=Psi+(data.y'-data.mu0)*(data.y-data.mu0'); nu_=nu+data.size;
    vecSigma_truesamp=zeros(N_prisamps,dim);
    for s=1:N_prisamps
        vecSigma_truesamp(s,:)=vech(iwishrnd(Psi_,nu_),'row');
    end
    % convert posterior samples and calculate estimates
    N_postsamps=length(resamp_idx);
    vecSigma_postsamp=zeros(N_postsamps,dim);
    sigma_sum=zeros(sqdim,1); Rho_sum=zeros(sqdim); Sigma_sum=zeros(sqdim);
%     for s=1:N_postsamps
%         [sigma_s,vecRho_s]=cnvrt_para(samp(s,1:sqdim)',samp(resamp_idx(s),sqdim+1:end),true);
%         sigma_sum=sigma_sum+sigma_s;
%         Rho_s=ivech(vecRho_s,'row'); Rho_s=Rho_s+tril(Rho_s,-1)';
%         Rho_sum=Rho_sum+Rho_s;
%         Sigma_s=sigma_s.*Rho_s.*sigma_s';
%         Sigma_sum=Sigma_sum+Sigma_s;
%         vecSigma_postsamp(s,:)=vech(Sigma_s,'row');
%     end
%     sigma_est=sigma_sum./N_postsamps; Rho_est=Rho_sum./N_postsamps; Sigma_est=Sigma_sum./N_postsamps;
    for s=1:N_postsamps
        [sigma_s,vecRho_s]=cnvrt_para(samp(s,1:sqdim)',samp(s,sqdim+1:end),true);
        sigma_sum=sigma_sum+wt(s).*sigma_s;
        Rho_s=ivech(vecRho_s,'row'); Rho_s=Rho_s+tril(Rho_s,-1)';
        Rho_sum=Rho_sum+wt(s).*Rho_s;
        Sigma_s=sigma_s.*Rho_s.*sigma_s';
        Sigma_sum=Sigma_sum+wt(s).*Sigma_s;
        vecSigma_postsamp(s,:)=vech(Sigma_s,'row');
    end
    sigma_est=sigma_sum./sum(wt); Rho_est=Rho_sum./sum(wt); Sigma_est=Sigma_sum./sum(wt);
    % plot posteriors
    for i=1:dim
        subplot(ha(i));
%         histogram(vecSigma_postsamp(:,i),'Normalization','pdf','displayname','$\Delta$-SphHMC'); hold on;
        [f,xi]=ksdensity(vecSigma_postsamp(:,i));
        plot(xi,f,'linestyle',lty{1},'linewidth',1.5,'displayname','$\Delta$-SphHMC'); hold on;
%         [f,xi]=ksdensity(vecSigma_truesamp(:,i));
        f=ksdensity(vecSigma_truesamp(:,i),xi);
        plot(xi,f,'linestyle',lty{2},'linewidth',1.5,'displayname','direct sampling'); hold on;
        set(gca,'fontsize',14);
        [I,J]=ind_vech2sub(sqdim,i,'row');
        if I==J
            xlim([0,1]);
            % add truth
            igampdf=@(X,A,B)gampdf(1./X,A,1./B)./(X.^2);
            f=igampdf(xi,(nu_-(sqdim-1))/2,Psi_(I,J)/2);
%             plot(xi,f,'linestyle',lty{3},'linewidth',1.5,'displayname','truth'); hold on;
        else
            xlim([-.5,.5]);
        end
        xlabel(['\sigma_{',num2str([I,J]),'}'],'fontsize',16);ylabel('Density','fontsize',16);
        h_lgd=legend('show','location','best');
        set(h_lgd,'fontsize',15,'interpreter','latex');
    end
end
% save estimates
% save(['./summary_bartlett/post_est',ifrewt,'.mat'],'sigma_est','Rho_est','Sigma_est');

fig.PaperPositionMode = 'auto';
print(fig,['./summary_bartlett/post_camp',ifrewt],'-dpng','-r0');