% This is to plot posterior density of correlations
% Prior is induced from squared-Dirichlet distribution

clear;
addpath('~/Documents/MATLAB/tight_subplot/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% plot for different parameters
alphas=[1,0.1,0.1];
L_alpha=length(alphas);
sqdim=3; dim=sqdim*(sqdim+1)/2;
% rho_sub=[1,2];
rho_sub=nchoosek(1:sqdim,2);
L_rho=size(rho_sub,1);
N_prisamps=1e6;

% setting
fig=figure(1); clf(fig);
% set(fig,'pos',[0 800 1000 350]);
% ha=tight_subplot(L_rho,L_alpha,[0.05,.06],[.14,.12],[.05,.05]);
set(fig,'pos',[0 800 1000 800]);
ha=tight_subplot(L_rho,L_alpha,[0.08,.05],[.06,.04],[.05,.04]);
lty={'-','--','-.',':'};

% plot priors
rho_idx=unique(rho_sub);
L=cell(sqdim,1); 
for i=1:L_alpha
    for k=rho_idx'
        % get samples from squared Dirichlet distribution
%         L{k}=sqdrchrnd(repmat(alpha(i),1,k),N_samps);
        L{k}=sqdrchrnd([repmat(alphas(i),1,k-1),1+9*(i==3)],N_prisamps);
    end
    for j=1:L_rho
        subplot(ha((j-1)*L_rho+i));
        d=min(rho_sub(j,:));
        samp_d=sum(L{rho_sub(j,1)}(:,1:d).*L{rho_sub(j,2)}(:,1:d),2);
        [f,xi]=ksdensity(samp_d(:),'support',[-1-1e-10,1+1e-10]);
        plot(xi,f,'linestyle',lty{1},'linewidth',1.5,'displayname','prior'); xlim([-1,1]);ylim([-.1,5]); hold on;
        set(gca,'fontsize',14);
        xlabel(['\rho_{',num2str(rho_sub(j,:)),'}'],'fontsize',14);ylabel('Density','fontsize',14);
        xlabh = get(gca,'XLabel');
        set(xlabh,'Position',get(xlabh,'Position') + [0 .1 0]);
%         title(['\alpha = ',num2str(alpha(i))],'fontsize',16);
        if i~=3
            title(['\alpha = (',repmat([num2str(alphas(i)),', '],1,sqdim-1),'1)'],'fontsize',16);
        else
            title(['\alpha = (',repmat([num2str(alphas(i)),', '],1,sqdim-1),'10)'],'fontsize',16);
        end
    end
end

% plot posteriors and MLEs
% load data
files = dir('./summary_sqdir');
nfiles = length(files) - 2;
sigma_est = zeros(sqdim,L_alpha);
Rho_est = zeros(sqdim,sqdim,L_alpha);
for i=1:L_alpha
    found=false;
    for j=1:nfiles
        if ~isempty(strfind(files(j+2).name,['_dim',num2str(dim),'_alpha',num2str(alphas(i)),'_']))
            if i==1 || contains(files(j+2).name,['_',num2str(10^(i-2)),'.mat'])
                load(strcat('./summary_sqdir/', files(j+2).name));
                fprintf('%s loaded.\n',files(j+2).name);
                found=true;
            end
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
% save('./summary_sqdir/post_est.mat','sigma_est','Rho_est');

fig.PaperPositionMode = 'auto';
print(fig,'./summary_sqdir/post_corr','-dpng','-r0');