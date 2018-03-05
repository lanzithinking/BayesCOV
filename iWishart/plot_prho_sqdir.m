% This is to plot prior density of correlations induced by
% squared-Dirichlet

clear;
addpath('~/Documents/MATLAB/tight_subplot/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% plot for different parameters
alpha=[1,0.5,0.1];
L_alpha=length(alpha);
rho_cmb=nchoosek(1:3,2);
L_rho=size(rho_cmb,1);
N_samps=1e6;

% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 350]);
ha=tight_subplot(1,L_alpha,[0,.07],[.12,.12],[.05,.05]);
lty={'-','-.','--',':','*','+'};

dim=max(rho_cmb(:));
rho_idx=unique(rho_cmb);
L=cell(dim,1); 
for i=1:L_alpha
    subplot(ha(i));
    for k=rho_idx'
        % get samples from squared Dirichlet distribution
%         L{k}=sqdrchrnd(repmat(alpha(i),1,k),N_samps);
        L{k}=sqdrchrnd([repmat(alpha(i),1,k-1),1],N_samps);
    end
    for j=1:L_rho
        d=min(rho_cmb(j,:));
        samp_d=sum(L{rho_cmb(j,1)}(:,1:d).*L{rho_cmb(j,2)}(:,1:d),2);
        [f,xi]=ksdensity(samp_d(:),'support',[-1-1e-10,1+1e-10]);
        plot(xi,f,'linestyle',lty{j},'linewidth',1.5);xlim([-1,1]);ylim([-.1,2]); hold on;
    end
    set(gca,'fontsize',14);
    xlabel('\rho','fontsize',14);ylabel('Density','fontsize',14);
    legend([repmat('\rho_{',L_rho,1),num2str(rho_cmb),repmat('}',L_rho,1)],'location','best');
%     title(['\alpha = ',num2str(alpha(i))],'fontsize',16);
    title(['\alpha = (',repmat([num2str(alpha(i)),', '],1,dim-1),'1)'],'fontsize',16);
end

fig.PaperPositionMode = 'auto';
print(fig,'./result/pri_corr_sqdir','-dpng','-r0');