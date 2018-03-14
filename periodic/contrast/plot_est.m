% This is to plot estimated mean and covariance functions

clear;
addpath('../../util/');
addpath('~/Documents/MATLAB/tight_subplot/');
addpath('~/Documents/MATLAB/boundedline/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% parameters for setting
M=20; N=100; sqdim=10; % dim=sqdim*(sqdim+1)/2;

% get posterior estimates
% set folder
folder = './summary/';
% load data
files = dir(folder);
nfiles = length(files) - 2;

found=false;
for j=1:nfiles
    if ~isempty(strfind(files(j+2).name,['_M',num2str(M),'_N',num2str(N),'_D',num2str(sqdim),'_'])) && contains(files(j+2).name,'_noreweight')
        load(strcat(folder, files(j+2).name));
        fprintf('%s loaded.\n',files(j+2).name);
        found=true;
    end
end
if found
%     N=length(t);
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
    samp_Sigma=samp_sigma.*samp_Rho.*reshape(samp_sigma,[Nsamp,N,1,sqdim]);
    Sigma_mean=squeeze(mean(samp_Sigma,1));
    Sigma_std=squeeze(std(samp_Sigma,0,1));
    % calculate the HPD
    mu_hpd=zeros(N,sqdim,2);
    Sigma_hpd=zeros(N,sqdim,sqdim,2);
    Rho_hpd=zeros(N,sqdim,sqdim,2);
    for n=1:N
        for i=1:sqdim
            mu_hpd(n,i,:)=calculate_hpd(samp_mu(:,n,i),0.95);
            for j=1:sqdim
                [Sigma_hpd(n,i,j,2),Sigma_hpd(n,i,j,1)]=calculate_hpd(samp_Sigma(:,n,i,j),0.95);
                [Rho_hpd(n,i,j,2),Rho_hpd(n,i,j,1)]=calculate_hpd(samp_Rho(:,n,i,j),0.95);
            end
        end
    end
end
% save estimates
filename=['periodic_sim_est_M',num2str(M),'_N',num2str(N),'_D',num2str(sqdim)];
save([folder,filename,'.mat'],...
    'seedNO','t','y','M_vecL','ker','a','b','m','V',...
    'alg_choice','Nsamp','NBurnIn','thin','stepsz','Nleap',...
    'mu_mean','mu_std','sigma_mean','sigma_std','Rho_mean','Rho_std','Sigma_mean','Sigma_std',...
    'mu_hpd','Rho_hpd','Sigma_hpd','acpt','time');

% load truth
load(['./periodic_multi_M',num2str(M),'_N',num2str(N),'_D',num2str(sqdim),'.mat'],'mu','Sigma','Rho');
mu=mu'; Sigma=permute(Sigma,[3,1,2]); Rho=permute(Rho,[3,1,2]);

% plot posterior estimates
% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 350]);
ha=tight_subplot(1,3,[0.1,.06],[.13,.08],[.055,.04]);
lty={'-','--','-.',':'};

% plot estimates
% plot mean
subplot(ha(1));
h1=plot(t,mu_mean-mu); hold on;
[l,p]=boundedline(t,mu_mean-mu,reshape(1.96.*mu_std,[size(mu_std,1),1,size(mu_std,2)]),'alpha'); hold on;
xlim([min(t),max(t)]); ylim([-3,3]);
set(gca,'box','on','fontsize',14);
xlabel('t','fontsize',16); ylabel('Error of Mean Estimate','fontsize',16);
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
xlabel('t','fontsize',16); ylabel('Error of Correlation Estimate','fontsize',16);
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