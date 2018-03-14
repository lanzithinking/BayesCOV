% This is to plot estimated mean and covariance functions from results of
% the latent factor process model by Fox and Dunson (2015). 

clear;
addpath('../../util/');
addpath('~/Documents/MATLAB/tight_subplot/');
addpath('~/Documents/MATLAB/boundedline/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% get posterior estimates
% set folder
folder = './summary/';
subfolder = [folder,'test/'];
% load settings
load([subfolder,'info4trial1.mat']);
% calculate estimates
[p,N,M]=size(y);
k = settings.k;
Niter = settings.Niter; Nburn = Niter/2;
saveEvery = settings.saveEvery;
trial = settings.trial;
% Build up some statistics of the samples

latent_mean = settings.latent_mean;
saveDir = [folder,settings.saveDir];

sampleEvery = settings.storeEvery;
Sigma_mean = zeros(N,p,p);
Sigma_std = zeros(N,p,p);
Sigma_hpd = zeros(N,p,p,2);
mu_mean = zeros(N,p);
mu_std = zeros(N,p);
mu_hpd = zeros(N,p,2);
Rho_mean = zeros(N,p,p);
Rho_std = zeros(N,p,p);
Rho_hpd = zeros(N,p,p,2);

cov_true = true_params.cov_true;
mu_true = true_params.mu;

for tt=1:N
    theta_zeta_tt = zeros(p,k,(Niter-Nburn)/sampleEvery);
    var_tt = zeros(p,p,(Niter-Nburn)/sampleEvery);
    mu_tt = zeros(p,(Niter-Nburn)/sampleEvery);
    m = 1;
    for nn=Nburn+1:sampleEvery:Niter
        n = nn+saveEvery-1;
        if rem(n,saveEvery)==0 & n<=Niter
            filename = [saveDir '/BNP_covreg_statsiter' num2str(n) 'trial' num2str(trial) '.mat'];
            load(filename)
            store_count = 1;
        end
        theta_zeta_tt(:,:,m) = Stats(store_count).theta*Stats(store_count).zeta(:,:,tt);
        var_tt(:,:,m) = Stats(store_count).theta*Stats(store_count).zeta(:,:,tt)*Stats(store_count).zeta(:,:,tt)'*Stats(store_count).theta'...
            + diag(1./Stats(store_count).invSig_vec);
        mu_tt(:,m) = Stats(store_count).theta*Stats(store_count).zeta(:,:,tt)*Stats(store_count).psi(:,tt);
        
        m = m + 1;
        store_count = store_count + 1;
    end
    sigma_tt = permute(var_tt,[3,1,2]); sigma_tt = sqrt(sigma_tt(:,1:p+1:end));
    Rho_tt = 1./reshape(sigma_tt',p,1,[]).*var_tt./reshape(sigma_tt',1,p,[]);
    
    Sigma_mean(tt,:,:) = mean(var_tt,3);
    Sigma_std(tt,:,:) = std(var_tt,0,3);
    
    Rho_mean(tt,:,:) = mean(Rho_tt,3);
    Rho_std(tt,:,:) = std(Rho_tt,0,3);
    
    mu_mean(tt,:) = mean(mu_tt,2);
    mu_std(tt,:) = std(mu_tt,0,2);
    
    for pp=1:p
        for jj=pp:p
            [Sigma_hpd(tt,pp,jj,2) Sigma_hpd(tt,pp,jj,1)] = calculate_hpd(var_tt(pp,jj,:),0.95);
            [Rho_hpd(tt,pp,jj,2) Rho_hpd(tt,pp,jj,1)] = calculate_hpd(Rho_tt(pp,jj,:),0.95);
        end
        if latent_mean
            [mu_hpd(tt,pp,2) mu_hpd(tt,pp,1)] = calculate_hpd(mu_tt(pp,:),0.95);
        end
    end
    
    if ~rem(tt,100)
        display(num2str(tt))
    end
    
end

% save estimates
filename=['periodic_sim_est_lfp_M',num2str(M),'_N',num2str(N),'_D',num2str(p)];
save([folder,filename,'.mat'],...
    'y','prior_params','settings','true_params',...
    'mu_mean','mu_std','mu_hpd','Sigma_mean','Sigma_std','Sigma_hpd','Rho_mean','Rho_std','Rho_hpd');

% load truth
load(['./periodic_multi_M',num2str(M),'_N',num2str(N),'_D',num2str(p),'.mat'],'mu','Sigma','Rho');
mu=mu'; Sigma=permute(Sigma,[3,1,2]); Rho=permute(Rho,[3,1,2]);

% plot posterior estimates
% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 350]);
ha=tight_subplot(1,3,[0.1,.06],[.13,.08],[.055,.04]);
lty={'-','--','-.',':'};

% plot estimates
sqdim=p; t=(1:N)'./N;
% plot mean
subplot(ha(1));
h1=plot(t,mu_mean-mu); hold on;
[l,p]=boundedline(t,mu_mean-mu,reshape(1.96.*mu_std,[size(mu_std,1),1,size(mu_std,2)]),'alpha'); hold on;
xlim([min(t),max(t)]); ylim([-3,3]);
set(gca,'box','on','fontsize',14);
xlabel('t','fontsize',16); ylabel('Error of Mean Estimaten','fontsize',16);
% set(gca,'colororderindex',1);
% plot(t,mu,'linestyle','--','linewidth',2); % add truth
if sqdim<=4
    lgd=legend(h1,[repmat('\mu_',sqdim,1),num2str((1:sqdim)')],'location','southwest');
    set(lgd,'orientation','horizontal','fontsize',14,'box','off');
end
% plot variance
subplot(ha(2));
idx=sub2ind([sqdim,sqdim],1:sqdim,1:sqdim);
h2=plot(t,Sigma_mean(:,idx)-Sigma(:,idx)); hold on;
[l,p]=boundedline(t,Sigma_mean(:,idx)-Sigma(:,idx),reshape(1.96.*Sigma_std(:,idx),[size(Sigma_std,1),1,length(idx)]),'alpha'); hold on;
xlim([min(t),max(t)]); ylim([-2,2]);
set(gca,'box','on','fontsize',14);
xlabel('t','fontsize',16); ylabel('Error of Variance Estimate','fontsize',16);
% set(gca,'colororderindex',1);
% plot(t,Sigma(:,idx),'linestyle','--','linewidth',2); % add truth
if sqdim<=4
    [I,J]=ind2sub([sqdim,sqdim],idx);
    lgd=legend(h2,[repmat('\sigma^2_{',length(I),1),num2str(I'),repmat('}',length(I),1)],'location','southeast');
    set(lgd,'orientation','horizontal','fontsize',14,'box','off');
end
% plot correlation
subplot(ha(3));
idx=find(logical(tril(ones(sqdim),-1)));
h3=plot(t,Rho_mean(:,idx)-Rho(:,idx)); hold on;
[l,p]=boundedline(t,Rho_mean(:,idx)-Rho(:,idx),reshape(1.96.*Rho_std(:,idx),[size(Rho_std,1),1,length(idx)]),'alpha'); hold on;
xlim([min(t),max(t)]); ylim([-1,1]);
set(gca,'box','on','fontsize',14);
xlabel('t','fontsize',16); ylabel('Estimated Correlation','fontsize',16);
% set(gca,'colororderindex',1);
% plot(t,Rho(:,idx),'linestyle','--','linewidth',2); % add truth
if sqdim<=4
    [I,J]=ind2sub([sqdim,sqdim],idx);
    lgd=legend(h3,[repmat('\rho_{',length(I),1),num2str(J),repmat(',',length(I),1),num2str(I),repmat('}',length(I),1)],'location','southeast');
    set(lgd,'orientation','horizontal','fontsize',14,'box','off');
end
% save plot
fig.PaperPositionMode = 'auto';
print(fig,[folder,filename],'-dpng','-r0');