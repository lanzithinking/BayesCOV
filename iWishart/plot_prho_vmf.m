% This is to plot prior density of correlations induced by von Mises-Fisher

clear;
addpath('~/Documents/MATLAB/tight_subplot/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% plot for different parameters
kappa=[1,10,100];
L_kappa=length(kappa);
rho_sub=nchoosek(1:3,2);
L_rho=size(rho_sub,1);
N_samps=1e6;

% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 350]);
ha=tight_subplot(1,L_kappa,[0,.07],[.12,.12],[.05,.05]);
lty={'-','-.','--',':','*','+'};

csvwrite('./tmp/mu.csv',eye(max(rho_sub(:))));
csvwrite('./tmp/rho_sub.csv',rho_sub);
for i=1:L_kappa
    subplot(ha(i));
    % generate samples by calling R script using movMF package
    csvwrite('./tmp/setting.csv',[N_samps,kappa(i)]);
    system('/usr/local/bin/R CMD BATCH gen_Rho_vmf.R ./tmp/output4debug.txt');
    Rho_i=csvread('./tmp/Rho.csv');
    for j=1:L_rho
        samp_j=Rho_i(:,j);
        [f,xi]=ksdensity(samp_j(:),'support',[-1-1e-10,1+1e-10]);
        plot(xi,f,'linestyle',lty{j},'linewidth',1.5);xlim([-1,1]);ylim([-.1,2]); hold on;
    end
    set(gca,'fontsize',14);
    xlabel('\rho','fontsize',14);ylabel('Density','fontsize',14);
    legend([repmat('\rho_{',L_rho,1),num2str(rho_sub),repmat('}',L_rho,1)],'location','best');
    title(['\kappa = ',num2str(kappa(i))],'fontsize',16);
end

fig.PaperPositionMode = 'auto';
print(fig,'./result/pri_corr_vmf','-dpng','-r0');