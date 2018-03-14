% This is to plot results contrasting our model and the latent factor
% process model by Fox and Dunson (2015).

clear;
addpath('../../util/');
addpath('~/Documents/MATLAB/tight_subplot/');
addpath('~/Documents/MATLAB/boundedline/');
addpath('~/Documents/MATLAB/supertitle/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% parameters for setting
M=10; N=200;
sqdim=2; %dim=sqdim*(sqdim+1)/2;

% get posterior estimates
% set folder
folder = {'./summary/','./test/'};
% plot posterior estimates
% setting
lty={'-','--','-.',':'};

% load truth
load(['./periodic_multi_M',num2str(M),'_N',num2str(N),'_D',num2str(sqdim),'.mat'],'mu','Sigma','Rho');
Sigma=permute(Sigma,[3,1,2]); Rho=permute(Rho,[3,1,2]);
% define filename
filename=cell(1,2);
filename{1}=['periodic_sim_est_M',num2str(M),'_N',num2str(N),'_D',num2str(sqdim)];
filename{2}=['periodic_sim_est_lfp_M',num2str(M),'_N',num2str(N),'_D',num2str(sqdim)];
for i=1:2
    % load estimates
    load([folder{i},[filename{i},'.mat']]);
    % plot posterior estimates
    fig=figure(i); clf(fig);
    set(fig,'pos',[0 800 800 350]);
    ha=tight_subplot(1,2,[0.1,.08],[.13,.11],[.07,.04]);
    % plot mean
    subplot(ha(1));
    h1=plot(t,mu_mean); hold on;
%     [l,p]=boundedline(t,mu_mean,reshape(1.96.*mu_std,[size(mu_std,1),1,size(mu_std,2)]),'alpha'); hold on;
    [l,p]=boundedline(t,mu_mean,permute(abs(mu_hpd(:,:,end:-1:1)-mu_mean),[1,3,2]),'alpha'); hold on;
    xlim([min(t),max(t)]); ylim([-1.5,2]);
    set(gca,'box','on','fontsize',16);
    xlabel('t','fontsize',20); ylabel('Estimated Mean','fontsize',20);
    set(gca,'colororderindex',1);
    plot(t,mu,'linestyle','--','linewidth',2); % add truth
    lgd=legend(h1,[repmat('\mu_',sqdim,1),num2str((1:sqdim)')],'location','southwest');
    set(lgd,'orientation','horizontal','fontsize',18,'box','off');
    % plot covariance
    subplot(ha(2));
    idx=[find(logical(eye(sqdim)));find(logical(tril(ones(sqdim),-1)))];
    h2=plot(t,Sigma_mean(:,idx)); hold on;
%     [l,p]=boundedline(t,Sigma_mean(:,idx),reshape(1.96.*Sigma_std(:,idx),[size(Sigma_std,1),1,length(idx)]),'alpha'); hold on;
    Sigma_hpd=permute(abs(Sigma_hpd(:,:,:,end:-1:1)-Sigma_mean),[1,4,2,3]);
    [l,p]=boundedline(t,Sigma_mean(:,idx),Sigma_hpd(:,:,idx),'alpha'); hold on;
    xlim([min(t),max(t)]); ylim([-.7,1.5]);
    set(gca,'box','on','fontsize',16);
    xlabel('t','fontsize',20); ylabel('Estimated Covariance','fontsize',20);
    set(gca,'colororderindex',1);
    plot(t,Sigma(:,idx),'linestyle','--','linewidth',2); % add truth
    [I,J]=ind2sub([sqdim,sqdim],idx);
    lgd=legend(h2,[repmat('\sigma_{',length(I),1),num2str(J),repmat(',',length(I),1),num2str(I),repmat('}',length(I),1)],'location','southwest');
    set(lgd,'orientation','horizontal','fontsize',18,'box','off');
    % suptitle(['M=',num2str(M),', N=',num2str(N)]);
    supertitle(['M=',num2str(M),', N=',num2str(N)],'fontsize',20);
    % save plot
    fig.PaperPositionMode = 'auto';
    print(fig,[folder{1},'contrast_',filename{i}],'-dpng','-r0');
end