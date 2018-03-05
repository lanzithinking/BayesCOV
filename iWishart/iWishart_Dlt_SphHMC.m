% simulation of inverse-Wishart posterior distribution using Dlt-SphHMC

clear;
addpath('../sampler/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% simulate data
sqdim=3; dim=sqdim*(sqdim+1)/2;
N=20;
mu0=zeros(sqdim,1);
% Psi=eye(sqdim); nu=sqdim+10;
% Sigma0=iwishrnd(Psi,nu);
a=1;
Sigma0=(eye(sqdim)+a)./11;
% generate data
y=mvnrnd(mu0',Sigma0,N);
% define data
data.size=N;
data.y=y; data.mu0=mu0; data.Sigma0=Sigma0;

% define priors
prior_tau.dim=sqdim;
prior_L.dim=dim;
%-- bartlett factor in inverse-wishart prior
Psi=eye(sqdim); nu=sqdim;
prior_tau.Psi=Psi; prior_tau.nu=nu;
prior_L.Psi=Psi; prior_L.nu=nu;
%-- (log-)normal prior
sigma2=0.01;
prior_tau.mean=zeros(sqdim,1);
prior_tau.var=sigma2.*ones(sqdim,1);
%-- squared-Dirichlet prior
alpha=.1;
prior_L.alpha=vech(alpha.*tril(ones(sqdim),-1)+10.*diag(ones(sqdim,1)),'row');
%-- von Mises_Fisher prior (polar direction)
kappa=100;
prior_L.kappa=kappa;
%-- Bingham prior (polar direction)
zeta=100;
prior_L.zeta=zeta;
% choices of priors
pri4tau_names={'bartlett','normal'};
pri4L_names={'bartlett','sqdir','vmf','bing'};
prior_tau.choice=pri4tau_names{2};
prior_L.choice=pri4L_names{2};
% print the info
fprintf('The prior for tau is set as %s; the prior for L is set as %s.\n',prior_tau.choice,prior_L.choice);

% sampling setting
stepsz=[1.2e-1,9e-4]; Nleap=[10,10];
alg={'standard','adaptive'};
alg_choice=alg{1};
reweight=false;

% allocation to save
Nsamp=1e6; burnrate=0.1; thin=10;
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
    vecL(idx_d)=rand(d,1);
    vecL(idx_d)=vecL(idx_d)./norm(vecL(idx_d));
end
% if contains(prior_L.choice,'bartlett')
%     Psi_=Psi+(data.y'-data.mu0)*(data.y-data.mu0'); nu_=nu+data.size;
%     nz_mat=1e-2.*randn(sqdim,1);
%     L=chol((Psi_(sqdim:-1:1,sqdim:-1:1)+nz_mat*nz_mat')./(nu_+sqdim+1),'lower');
%     L=L./sqrt(sum(L.^2,2));
%     vecL=vech(L,'row');
% end
% define Dlt-SphHMC sampler
dlt_sphHMC4L=Delta_SphHMC(vecL,@(vecL,adj_metric)geom_L(vecL,prior_L,data,tau,adj_metric),stepsz(2),Nleap(2),alg_choice,reweight);

fprintf('Running %s Delta-Spherical HMC with %s reweight...\n',alg_choice,~reweight*'no');
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
    % update potential and its gradient for tau
    [U_tau,dU_tau]=geom_tau(tau,prior_tau,data,vecL);
    [tau,U_tau,dU_tau,acpt_tau]=HMC(tau,U_tau,dU_tau,@(tau)geom_tau(tau,prior_tau,data,vecL),stepsz(1),Nleap(1));
    
    % Use Delta-Spherical HMC to sample L
    % update geometry function for L
    dlt_sphHMC4L.geom=@(vecL)geom_L(vecL,prior_L,data,tau,~dlt_sphHMC4L.reweight);
    % update potential and its gradient for L
    [dlt_sphHMC4L.U,dlt_sphHMC4L.dU]=dlt_sphHMC4L.geom(vecL);
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
    case 'bing'
        sufile=strcat(['_zeta',num2str(prior_L.zeta)]);
    otherwise
        error('Prior not defined!');
end
if reweight
    ifrewt='_reweight';
else
    ifrewt='_noreweight';
end
f_name=['iWishart_dlt_SphHMC_dim',num2str(dim),sufile,'_',time_lbl,ifrewt];
if exist('result','dir')~=7
    mkdir('result');
end
save(['result/',f_name,'.mat'],'data','prior_tau','prior_L','alg_choice','Nsamp','NBurnIn','stepsz','Nleap','samp','reweight','wt','resamp_idx','acpt','time');
% summarize
fprintf('\nIt takes %.2f seconds to collect %e samples after thinning %d.\n', time,Nsamp,thin);
fprintf('\nThe final acceptance rate is: %.4f \t %.4f \n\n',acpt);
% efficiency measurement
addpath('./result/');
CalculateStatistics(f_name,'./result/');