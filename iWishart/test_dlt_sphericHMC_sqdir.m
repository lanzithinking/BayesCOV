% This is to test Delta-Spherical HMC

clear;
addpath('~/Documents/MATLAB/tight_subplot/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% plot for different parameters
alphas=0.4;
L_alpha=length(alphas);
sqdim=3; dim=sqdim*(sqdim+1)/2;
% rho_cmb=[1,2];
rho_cmb=nchoosek(1:sqdim,2);
L_rho=size(rho_cmb,1);
N_prisamps=1e6;
reweight=false;
if reweight
    ifrewt='_reweight';
else
    ifrewt='_noreweight';
end

% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 350]);
ha=tight_subplot(L_alpha,L_rho,[0.05,.06],[.14,.12],[.05,.05]);
% set(fig,'pos',[0 800 1000 800]);
% ha=tight_subplot(L_alpha,L_rho,[0.08,.05],[.06,.05],[.05,.05]);
lty={'-','--','-.',':'};

% plot samples from direct sampling
rho_idx=unique(rho_cmb);
L=cell(sqdim,1);
for i=1:L_alpha
    for k=rho_idx'
        % get samples from squared Dirichlet distribution
%         L{k}=sqdrchrnd(repmat(alpha(i),1,k),N_samps);
        L{k}=sqdrchrnd([repmat(alphas(i),1,k-1),1],N_prisamps);
    end
    for j=1:L_rho
        subplot(ha((i-1)*L_alpha+j));
        d=min(rho_cmb(j,:));
        samp_d=sum(L{rho_cmb(j,1)}(:,1:d).*L{rho_cmb(j,2)}(:,1:d),2);
        [f,xi]=ksdensity(samp_d(:),'support',[-1-1e-10,1+1e-10]);
        plot(xi,f,'linestyle',lty{mod(i-1,length(lty))+1},'linewidth',1.5,'displayname','truth'); xlim([-1,1]);ylim([-.1,2]); hold on;
        set(gca,'fontsize',14);
        xlabel(['\rho_{',num2str(rho_cmb(j,:)),'}'],'fontsize',14);ylabel('Density','fontsize',14);
%         title(['\alpha = ',num2str(alpha(i))],'fontsize',16);
        title(['\alpha = (',repmat([num2str(alphas(i)),', '],1,sqdim-1),'1)'],'fontsize',16);
    end
end

% plot samples from delta-spherical HMC
% load data
files = dir('./summary_sqdir');
nfiles = length(files) - 2;
for i=1:L_alpha
    found=false;
    for j=1:nfiles
        if contains(files(j+2).name, 'test_') && contains(files(j+2).name, 'ifrewt') && ~isempty(strfind(files(j+2).name,['_dim',num2str(dim),'_alpha',num2str(alphas(i)),'_'])) 
            load(strcat('./summary_sqdir/', files(j+2).name));
            found=true;
        end
    end
    if ~found
        addpath('../sampler/');
        % simulate data
        sqdim=3; dim=sqdim*(sqdim+1)/2;
        N=20;
        mu0=zeros(sqdim,1);
        a=1;
        Sigma0=(eye(sqdim)+a)./11;
        % generate data
        y=zeros(N,sqdim); % force zero data to sample from the prior
        % define data
        data.size=N;
        data.y=y; data.mu0=mu0; data.Sigma0=Sigma0;

        % define priors
        prior_tau.dim=sqdim;
        prior_L.dim=dim;
        %-- bartlett factor in inverse-wishart prior
        Psi=eye(sqdim); nu=sqdim+10;
        prior_tau.Psi=Psi; prior_tau.nu=nu;
        prior_L.Psi=Psi; prior_L.nu=nu;
        %-- (log-)normal prior
        sigma2=0.01;
        prior_tau.mean=zeros(sqdim,1);
        prior_tau.var=sigma2.*ones(sqdim,1);
        %-- squared-Dirichlet prior
%         alphas=0.4;
        prior_L.alpha=vech(alphas(i).*tril(ones(sqdim),-1)+diag(ones(sqdim,1)),'row');
        %-- von Mises_Fisher prior
        kappa=100;
        prior_L.kappa=kappa;
        % choices of priors
        pri4tau_names={'bartlett','normal'};
        pri4L_names={'bartlett','sqdir','vmf'};
        prior_tau.choice=pri4tau_names{2};
        prior_L.choice=pri4L_names{2};
        % print the info
        fprintf('The prior for tau is set as %s; the prior for L is set as %s.\n',prior_tau.choice,prior_L.choice);

        % sampling setting
        stepsz=[1e-1,1e-1]; Nleap=[10,10];
        alg={'standard','robust'};
        alg_choice=alg{1};
%         reweight=false;

        % allocation to save
        Nsamp=1e6; burnrate=0.1; thin=1;
        Niter=Nsamp*thin; NBurnIn=floor(Niter*burnrate); Niter=Niter+ NBurnIn;
        samp=zeros(Nsamp,sqdim+dim);
        engy=zeros(Niter,2);
        acpt=zeros(1,2); % overall acceptance
        accp=zeros(1,2); % online acceptance

        % initializatioin
        tau=1e-2.*randn(sqdim,1);
        vecL=zeros(dim,1);
        for d=1:sqdim
            idx_d=1+d*(d-1)/2:d*(d+1)/2;
            vecL(idx_d)=randn(d,1);
            vecL(idx_d)=vecL(idx_d)./norm(vecL(idx_d));
        end
        [U_tau,dU_tau]=geom_tau(tau,prior_tau,data,vecL);
        % define Dlt-SphHMC sampler
        dlt_sphHMC4L=Delta_SphHMC(vecL,@(vecL,adj_metric)geom_L(vecL,prior_L,data,tau,adj_metric),stepsz(2),Nleap(2),alg_choice,reweight);

        fprintf('Running %s Delta-Spherical HMC with %sreweight...\n',alg_choice,~reweight*'no');
        prog=0.05:0.05:1;
        tic;
        for iter=1:Niter

            % display sampleing progress and online acceptance rate
            if ismember(iter,floor(Niter.*prog))
                fprintf('%.0f%% iterations completed.\n',100*iter/Niter);
                fprintf('Online acceptance rates: %.0f%%, %.0f%%\n', accp.*(100/floor(prog(1)*Niter)));
                accp=zeros(1,2);
            end

            % Use HMC to sample tau
            [tau,U_tau,dU_tau,acpt_tau]=HMC(tau,U_tau,dU_tau,@(tau)geom_tau(tau,prior_tau,data,vecL),stepsz(1),Nleap(1));
%             acpt_tau=true;

            % Use Delta-Spherical HMC to sample L
            % update geometry function for L
            dlt_sphHMC4L.geom=@(vecL)geom_L(vecL,prior_L,data,tau,~dlt_sphHMC4L.reweight);
            [dlt_sphHMC4L,acpt_L]=dlt_sphHMC4L.std_Dlt_SphHMC(true);
            vecL=dlt_sphHMC4L.q;

            % accptances
            accp=accp+[acpt_tau,acpt_L];

            % burn-in complete
            if(iter==NBurnIn)
                fprintf('Burn in completed!\n');
            end

            engy(iter,:)=[U_tau,dlt_sphHMC4L.U];
            % save samples after burn-in
            if(iter>NBurnIn)
                if mod(iter-NBurnIn,thin) == 0
                    samp(ceil((iter-NBurnIn)/thin),:)=[tau',vecL'];
                end
                acpt=acpt+[acpt_tau,acpt_L];
            end

        end

        % count time
        time=toc;
        % save results
        acpt=acpt/(Niter-NBurnIn);
        % resample if set to reweight
        if reweight
            wt=sum(log(abs(samp(:,sqdim+(1:sqdim).*((1:sqdim)+1)./2))),2);
            wt=exp(wt-median(wt));
            resamp_idx=datasample(1:Nsamp,Nsamp,'replace',true,'weights',wt);
        else
            wt=ones(Nsamp,1);
            resamp_idx=1:Nsamp;
        end
        % save
        time_lbl=regexprep(num2str(fix(clock)),'    ','_');
        switch lower(prior_L.choice)
            case 'bartlett'
                sufile=strcat(['_nu',num2str(prior_L.nu)]);
            case 'sqdir'
                sufile=strcat(['_alpha',num2str(prior_L.alpha(2))]);
            case 'vmf'
                sufile=strcat(['_kappa',num2str(prior_L.kappa)]);
            otherwise
                error('Prior not defined!');
        end
        f_name=['test_iWishart_dlt_SphHMC_dim',num2str(dim),sufile,'_',time_lbl,ifrewt];
        save(['summary_sqdir/',f_name,'.mat'],'data','prior_tau','prior_L','alg_choice','Nsamp','NBurnIn','stepsz','Nleap','samp','wt','resamp_idx','acpt','time');
        % summarize
        fprintf('\nIt takes %.2f seconds to collect %e samples after thinning %d.\n', time,Nsamp,thin);
        fprintf('\nThe final acceptance rate is: %.4f \t %.4f \n\n',acpt);
    end
    % convert samples and calculate estimates
    N_postsamp=length(resamp_idx);
    rho=zeros(L_rho,N_postsamp);
    for s=1:N_postsamp
        [sigma,vecRho]=cnvrt_para(samp(s,1:sqdim)',samp(resamp_idx(s),sqdim+1:end));
        rho(:,s)=vecRho(sub2ind_vech(sqdim,rho_cmb(:,2),rho_cmb(:,1),'row'));
    end
    for j=1:L_rho
        subplot(ha((i-1)*L_alpha+j));
        [f,xi]=ksdensity(rho(j,:),'support',[-1-1e-10,1+1e-10]);
        plot(xi,f,'linestyle',lty{mod(i-1,length(lty))+2},'linewidth',1.5,'displayname','estimate'); xlim([-1,1]);ylim([-.1,2]); hold on;
        legend('show','location','best');
    end
end

fig.PaperPositionMode = 'auto';
print(fig,'./summary_sqdir/test_dlt_sphericalHMC','-dpng','-r0');