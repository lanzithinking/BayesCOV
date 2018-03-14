% This is to plot estimated mean and covariance functions

clear;
addpath('../util/');
addpath('~/Documents/MATLAB/tight_subplot/');
addpath('~/Documents/MATLAB/boundedline/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% parameters for setting
sqdim=2; dim=sqdim*(sqdim+1)/2;

% get posterior estimates using Dlt-Spherical HMC
% load data
files = dir('./result');
nfiles = length(files) - 2;

found=false;
for j=1:nfiles
    if ~isempty(strfind(files(j+2).name,['_dim',num2str(dim),'_'])) && contains(files(j+2).name,'_noreweight')
        load(strcat('./result/', files(j+2).name));
        fprintf('%s loaded.\n',files(j+2).name);
        found=true;
    end
end
if found
    N=length(t);
    % calculate estimates with posterior samples
    mu_mean=squeeze(mean(samp_mu,1));
    mu_std=squeeze(std(samp_mu,0,1));
    samp_sigma=exp(samp_tau);
    sigma_mean=squeeze(mean(samp_sigma,1));
    sigma_std=squeeze(std(samp_sigma,0,1));
    samp_Rho=zeros(Nsamp,N,sqdim,sqdim);
    for i=1:sqdim
        for j=1:sqdim
            ind_i=1+i*(i-1)/2:i*(i+1)/2;
            ind_j=1+j*(j-1)/2:j*(j+1)/2;
            min_ij=min([i,j]);
            samp_Rho(:,:,i,j)=sum(samp_vecL(:,:,ind_i(1:min_ij)).*samp_vecL(:,:,ind_j(1:min_ij)),3);
        end
    end
    Rho_mean=squeeze(mean(samp_Rho,1));
    Rho_std=squeeze(std(samp_Rho,0,1));
    Sigma_mean=squeeze(mean(samp_sigma.*samp_Rho.*reshape(samp_sigma,[Nsamp,N,1,sqdim]),1));
    Sigma_std=squeeze(std(samp_sigma.*samp_Rho.*reshape(samp_sigma,[Nsamp,N,1,sqdim]),0,1));
end
% save estimates
save('./result/periodic_sim_est.mat',...
    'seedNO','t','y','M_vecL','ker','a','b','m','V',...
    'alg_choice','Nsamp','NBurnIn','thin','stepsz','Nleap',...
    'mu_mean','mu_std','sigma_mean','sigma_std','Rho_mean','Rho_std','Sigma_mean','Sigma_std','acpt','time');

% load truth
load('./periodic_data.mat','mu','Sigma','Rho');
Sigma=permute(Sigma,[3,1,2]); Rho=permute(Rho,[3,1,2]);

% plot posterior estimates
% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 350]);
nrow=floor(sqrt(dim)); ncol=ceil(dim/nrow);
ha=tight_subplot(1,3,[0.1,.06],[.13,.08],[.055,.04]);
lty={'-','--','-.',':'};

% plot estimates
% plot mean
subplot(ha(1));
h1=plot(t,mu_mean); hold on;
[l,p]=boundedline(t,mu_mean,reshape(1.96.*mu_std,[size(mu_std,1),1,size(mu_std,2)]),'alpha'); hold on;
xlim([min(t),max(t)]);
set(gca,'box','on','fontsize',14);
xlabel('t','fontsize',16); ylabel('Estimated Mean','fontsize',16);
set(gca,'colororderindex',1);
plot(t,mu,'linestyle','--','linewidth',2); % add truth
lgd=legend(h1,[repmat('\mu_',sqdim,1),num2str((1:sqdim)')],'location','southwest');
set(lgd,'fontsize',14);
% plot covariance
subplot(ha(2));
h2=plot(t,Sigma_mean(:,[1,4,2])); hold on;
[l,p]=boundedline(t,Sigma_mean(:,[1,4,2]),reshape(1.96.*Sigma_std(:,[1,4,2]),[size(Sigma_std,1),1,3]),'alpha'); hold on;
xlim([min(t),max(t)]); % ylim([-0.5,2]);
set(gca,'box','on','fontsize',14);
xlabel('t','fontsize',16); ylabel('Estimated Variance/Covariance','fontsize',16);
set(gca,'colororderindex',1);
plot(t,Sigma(:,[1,4,2]),'linestyle','--','linewidth',2); % add truth
[I,J]=ind_vech2sub(sqdim,[1,3,2],'row');
lgd=legend(h2,[repmat('\sigma_{',length(I),1),num2str(J),repmat(',',length(I),1),num2str(I),repmat('}',length(I),1)],'location','southeast');
set(lgd,'orientation','horizontal','fontsize',14);
% plot correlation
subplot(ha(3));
h3=plot(t,Rho_mean(:,2)); hold on;
[l,p]=boundedline(t,Rho_mean(:,2),reshape(1.96.*Rho_std(:,2),[size(Rho_std,1),1,1]),'alpha'); hold on;
xlim([min(t),max(t)]);
set(gca,'box','on','fontsize',14);
xlabel('t','fontsize',16); ylabel('Estimated Correlation','fontsize',16);
set(gca,'colororderindex',1);
plot(t,Rho(:,2),'linestyle','--','linewidth',2); % add truth
[I,J]=ind_vech2sub(sqdim,2,'row');
lgd=legend(h3,[repmat('\rho_{',length(I),1),num2str(J),repmat(',',length(I),1),num2str(I),repmat('}',length(I),1)],'location','southeast');
set(lgd,'fontsize',14);
% save plot
fig.PaperPositionMode = 'auto';
print(fig,'./result/periodic_sim_est','-dpng','-r0');