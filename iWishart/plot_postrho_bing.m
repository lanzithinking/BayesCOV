% This is to plot posterior density of correlations
% Prior is induced from Bingham distribution

clear;
addpath('~/Documents/MATLAB/tight_subplot/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% plot for different parameters
zeta=[1,10,100];
L_zeta=length(zeta);
sqdim=3; dim=sqdim*(sqdim+1)/2;
% rho_sub=[1,2];
rho_sub=nchoosek(1:sqdim,2);
L_rho=size(rho_sub,1);
N_prisamps=1e6;

% setting
fig=figure(1); clf(fig);
% set(fig,'pos',[0 800 1000 350]);
% ha=tight_subplot(L_rho,L_zeta,[0.05,.06],[.14,.12],[.05,.05]);
set(fig,'pos',[0 800 1000 800]);
ha=tight_subplot(L_rho,L_zeta,[0.08,.05],[.06,.04],[.05,.04]);
lty={'-','--','-.',':'};

% plot priors
csvwrite('./tmp/mu.csv',eye(max(rho_sub(:))));
csvwrite('./tmp/rho_sub.csv',rho_sub);
for i=1:L_zeta
    % generate samples by calling R script using Directional package
    csvwrite('./tmp/setting.csv',[N_prisamps,zeta(i)]);
    system('/usr/local/bin/R CMD BATCH gen_Rho_bing.R ./tmp/output4debug.txt');
    Rho_i=csvread('./tmp/Rho.csv');
    for j=1:L_rho
        subplot(ha((j-1)*L_rho+i));
        samp_j=Rho_i(:,j);
        [f,xi]=ksdensity(samp_j(:),'support',[-1-1e-10,1+1e-10]);
        plot(xi,f,'linestyle',lty{1},'linewidth',1.5,'displayname','prior'); xlim([-1,1]);ylim([-.1,5]); hold on;
        set(gca,'fontsize',14);
        xlabel(['\rho_{',num2str(rho_sub(j,:)),'}'],'fontsize',14);ylabel('Density','fontsize',14);
        xlabh = get(gca,'XLabel');
        set(xlabh,'Position',get(xlabh,'Position') + [0 .1 0]);
        title(['\zeta = ',num2str(zeta(i))],'fontsize',16);
    end
end

% plot posteriors and MLEs
% load data
files = dir('./summary_bing');
nfiles = length(files) - 2;
sigma_est = zeros(sqdim,L_zeta);
Rho_est = zeros(sqdim,sqdim,L_zeta);
for i=1:L_zeta
    found=false;
    for j=1:nfiles
        if ~isempty(strfind(files(j+2).name,['_dim',num2str(dim),'_zeta',num2str(zeta(i)),'_']))
            load(strcat('./summary_bing/', files(j+2).name));
            fprintf('%s loaded.\n',files(j+2).name);
            found=true;
        end
    end
    if found
        % convert posterior samples and calculate estimates
        N_postsamp=length(resamp_idx);
        rho=zeros(L_rho,N_postsamp);
        sigma_sum=zeros(sqdim,1); vecRho_sum=zeros(dim,1);
        idx_vech=sub2ind_vech(sqdim,rho_sub(:,2),rho_sub(:,1),'row');
        for s=1:N_postsamp
            [sigma,vecRho]=cnvrt_para(samp(s,1:sqdim)',samp(resamp_idx(s),sqdim+1:end));
            rho(:,s)=vecRho(idx_vech);
            sigma_sum=sigma_sum+sigma;
            vecRho_sum=vecRho_sum+vecRho;
        end
        sigma_est(:,i)=sigma_sum./N_postsamp;
        Rho_est(:,:,i)=ivech(vecRho_sum./N_postsamp,'row');
        % obtain sample correlation coefficients
        R=corr(data.y);
        r=R(sub2ind([sqdim,sqdim],rho_sub(:,1),rho_sub(:,2)));
        for j=1:L_rho
            subplot(ha((j-1)*L_rho+i));
            [f,xi]=ksdensity(rho(j,:),'support',[-1-1e-10,1+1e-10]);
            plot(xi,f,'linestyle',lty{2},'linewidth',1.5,'displayname','posterior'); xlim([-1,1]);ylim([-.1,5]); hold on;
            plot([r(j),r(j)],ylim,lty{3},'linewidth',1.5,'displayname','MLE'); hold on;
            legend('show','location','northwest');
        end
    end
end
% save estimates
save('./summary_bing/post_est.mat','sigma_est','Rho_est');

fig.PaperPositionMode = 'auto';
print(fig,'./summary_bing/post_corr','-dpng','-r0');