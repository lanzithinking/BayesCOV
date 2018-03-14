% This is to plot simulated data

clear;
addpath('../../util/');
addpath('~/Documents/MATLAB/tight_subplot/');
addpath('~/Documents/MATLAB/columnlegend/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% parameters for setting
M=100; N=200;
sqdim=2; %dim=sqdim*(sqdim+1)/2;

% load truth
load(['./periodic_multi_M',num2str(M),'_N',num2str(N),'_D',num2str(sqdim),'.mat']);
mu=mu'; Sigma=permute(Sigma,[3,1,2]); Rho=permute(Rho,[3,1,2]);

% plot posterior estimates
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 800 350]);
ha=tight_subplot(1,2,[0.1,.08],[.13,.11],[.07,.05]);
lty={'-','--','-.',':'};

% plot y_1 and mean
subplot(ha(1));
plot(t,squeeze(y(1,:,:)),'*','markersize',4.5); hold on;
set(gca,'colororderindex',1);
plot(t,mu,'linewidth',2);
xlim([min(t),max(t)]); ylim([-3,3]);
set(gca,'fontsize',14);
xlabel('t','fontsize',16); ylabel('y','rot',0,'fontsize',16);
lgd_str=cell(2*D,1);
lgd_str(1:D)=cellstr([repmat('y_',D,1),num2str((1:D)')]);
lgd_str(1+D:end)=cellstr([repmat('\mu_',D,1),num2str((1:D)')]);
lgd=columnlegend(D,lgd_str,'location','south');
set(lgd,'orientation','horizontal','fontsize',15,'box','off');
% plot covariance
subplot(ha(2));
idx=[find(logical(eye(sqdim)));find(logical(tril(ones(sqdim),-1)))];
h2=plot(t,Sigma(:,idx)); hold on;
xlim([min(t),max(t)]); ylim([-.5,1]);
set(gca,'box','on','fontsize',14);
xlabel('t','fontsize',16); ylabel('Covariance','fontsize',16);
[I,J]=ind2sub([sqdim,sqdim],idx);
lgd=legend(h2,[repmat('\sigma_{',length(I),1),num2str(J),repmat(',',length(I),1),num2str(I),repmat('}',length(I),1)],'location','southeast');
set(lgd,'orientation','horizontal','fontsize',15,'box','off');
% save plot
fig.PaperPositionMode = 'auto';
print(fig,'periodic_sim_data','-dpng','-r0');
