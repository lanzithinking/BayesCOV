% This is to compare selected mean and covariance functions between truth
% and estimates

clear;
addpath('../util/');
addpath('~/Documents/MATLAB/tight_subplot/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% parameters for setting
N=100;
sqdim=10; dim=sqdim*(sqdim+1)/2;

% load truth
if exist('./LFp_data.mat','file')
    % load data
    load('./LFp_data.mat','mu','Sigma');
    mu_truth=mu';
    Sigma_truth=permute(Sigma,[3,1,2]);
end

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
    % calculate estimates with posterior samples
    mu_est=squeeze(mean(samp_mu,1));
    mu_std=squeeze(std(samp_mu));
    sigma_est=squeeze(mean(exp(samp_tau),1));
    sigma_std=squeeze(std(exp(samp_tau)));
    samp_Rho=zeros(Nsamp,N,sqdim,sqdim);
    for i=1:sqdim
        for j=1:sqdim
            ind_i=1+i*(i-1)/2:i*(i+1)/2;
            ind_j=1+j*(j-1)/2:j*(j+1)/2;
            min_ij=min([i,j]);
            samp_Rho(:,:,i,j)=sum(samp_vecL(:,:,ind_i(1:min_ij)).*samp_vecL(:,:,ind_j(1:min_ij)),3);
        end
    end
    Rho_est=squeeze(mean(samp_Rho,1));
    Rho_std=squeeze(std(samp_Rho));
    Sigma_est=squeeze(mean(exp(samp_tau).*samp_Rho.*reshape(exp(samp_tau),[Nsamp,N,1,sqdim]),1));
    Sigma_std=squeeze(std(exp(samp_tau).*samp_Rho.*reshape(exp(samp_tau),[Nsamp,N,1,sqdim])));
end
% save estimates
save('./result/LFp_estd.mat','mu_est','sigma_est','Rho_est','Sigma_est',...
                             'mu_std','sigma_std','Rho_std','Sigma_std');

% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 600]);
nrow=floor(sqrt(dim)); ncol=ceil(dim/nrow);
ha=tight_subplot(3,3,[0.08,.06],[.08,.05],[.06,.04]);
lty={'-','--','-.',':'};

% selected means
sel_mu=[1,floor(sqdim/2),sqdim];
% selected covariance
sel_cov=5:sqdim:dim;

% plot posterior estimates
% plot mean
for i=1:length(sel_mu)
    subplot(ha(i));
    plot(t,mu_truth(:,sel_mu(i)),'color','red','linestyle',lty{1},'linewidth',2); hold on;
    plot(t,mu_est(:,sel_mu(i)),'color','blue','linestyle',lty{1},'linewidth',1.5); hold on;
    plot(t,[mu_est(:,sel_mu(i))-1.96*mu_std(:,sel_mu(i)),mu_est(:,sel_mu(i))+1.96*mu_std(:,sel_mu(i))],...
        'color','black','linestyle',lty{2},'linewidth',1.5); hold on;
    set(gca,'fontsize',14);
    xlabel('x','fontsize',16); ylabel(['$\mu_{',num2str(sel_mu(i)),'}$'],'interpreter','latex','fontsize',16);
end
% plot covariance
[I,J]=ind_vech2sub(sqdim,sel_cov,'row');
for i=1:length(sel_cov)
    subplot(ha(length(sel_mu)+i));
    plot(t,Sigma_truth(:,I(i),J(i)),'color','red','linestyle',lty{1},'linewidth',2); hold on;
    plot(t,Sigma_est(:,I(i),J(i)),'color','blue','linestyle',lty{1},'linewidth',1.5); hold on;
    plot(t,[Sigma_est(:,I(i),J(i))-1.96*Sigma_std(:,I(i),J(i)),Sigma_est(:,I(i),J(i))+1.96*Sigma_std(:,I(i),J(i))],...
        'color','black','linestyle',lty{2},'linewidth',1.5); hold on;
    set(gca,'fontsize',14);
    xlabel('x','fontsize',16); ylabel(['$\Sigma_{',num2str(I(i)),',',num2str(J(i)),'}$'],'interpreter','latex','fontsize',16);
end

% save plots
fig.PaperPositionMode = 'auto';
print(fig,'./result/LFp_sim_estd','-dpng','-r0');