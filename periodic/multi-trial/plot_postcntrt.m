% This is to plot posterior contraction phenomena

clear;
addpath('../../util/');
addpath('~/Documents/MATLAB/tight_subplot/');
addpath('~/Documents/MATLAB/boundedline/');
addpath('~/Documents/MATLAB/supertitle/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% parameters for setting
M=[10,100]; N=[20,200];
n_M=length(M); n_N=length(N);
sqdim=2; %dim=sqdim*(sqdim+1)/2;

% get posterior estimates
% set folder
folder = './summary/';
% plot posterior estimates
% setting
lty={'-','--','-.',':'};
for i=1:n_M
    for j=1:n_N
        % load estimates
        load([folder,['periodic_sim_est_M',num2str(M(i)),'_N',num2str(N(j)),'_D',num2str(sqdim),'.mat']])
        % load truth
        load(['./periodic_multi_M',num2str(M(i)),'_N',num2str(N(j)),'_D',num2str(sqdim),'.mat'],'mu','Sigma','Rho');
        Sigma=permute(Sigma,[3,1,2]); Rho=permute(Rho,[3,1,2]);
        % plot posterior estimates
        fig=figure((i-1)*n_M+j); clf(fig);
        set(fig,'pos',[0 800 800 350]);
        ha=tight_subplot(1,2,[0.1,.08],[.13,.11],[.07,.04]);
        % plot mean
        subplot(ha(1));
        h1=plot(t,mu_mean); hold on;
        [l,p]=boundedline(t,mu_mean,reshape(1.96.*mu_std,[size(mu_std,1),1,size(mu_std,2)]),'alpha'); hold on;
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
        [l,p]=boundedline(t,Sigma_mean(:,idx),reshape(1.96.*Sigma_std(:,idx),[size(Sigma_std,1),1,length(idx)]),'alpha'); hold on;
        xlim([min(t),max(t)]); ylim([-.5,1.5]);
        set(gca,'box','on','fontsize',16);
        xlabel('t','fontsize',20); ylabel('Estimated Covariance','fontsize',20);
        set(gca,'colororderindex',1);
        plot(t,Sigma(:,idx),'linestyle','--','linewidth',2); % add truth
        [I,J]=ind2sub([sqdim,sqdim],idx);
        lgd=legend(h2,[repmat('\sigma_{',length(I),1),num2str(J),repmat(',',length(I),1),num2str(I),repmat('}',length(I),1)],'location','southwest');
        set(lgd,'orientation','horizontal','fontsize',18,'box','off');
%         suptitle(['M=',num2str(M(i)),', N=',num2str(N(j))]);
        supertitle(['M=',num2str(M(i)),', N=',num2str(N(j))],'fontsize',20);
        % save plot
        fig.PaperPositionMode = 'auto';
        print(fig,[folder,['periodic_postcntrt_M',num2str(M(i)),'_N',num2str(N(j))]],'-dpng','-r0');
    end
end