% this is to test the calculation of gradients

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
Psi=eye(sqdim); nu=sqdim+10;
prior_tau.Psi=Psi; prior_tau.nu=nu;
prior_L.Psi=Psi; prior_L.nu=nu;
%-- (log-)normal prior
sigma2=0.01;
prior_tau.mean=zeros(sqdim,1);
prior_tau.var=sigma2.*ones(sqdim,1);
%-- squared-Dirichlet prior
alpha=0.1;
prior_L.alpha=vech(alpha.*tril(ones(sqdim),-1)+diag(ones(sqdim,1)),'row');
%-- von Mises_Fisher prior (polar direction)
kappa=100;
prior_L.kappa=kappa;
%-- Bingham prior (polar direction)
zeta=100;
prior_L.zeta=zeta;
% choices of priors
pri4tau_names={'bartlett','normal'};
pri4L_names={'bartlett','sqdir','vmf','bing'};
prior_tau.choice=pri4tau_names{1};
prior_L.choice=pri4L_names{1};
% print the info
fprintf('The prior for tau is set as %s; the prior for L is set as %s.\n',prior_tau.choice,prior_L.choice);

% calculation setting
reweight=false;

% initializatioin
tau=1e-2.*randn(sqdim,1);
vecL=zeros(dim,1);
for d=1:sqdim
    idx_d=1+d*(d-1)/2:d*(d+1)/2;
    vecL(idx_d)=randn(d,1);
    vecL(idx_d)=vecL(idx_d)./norm(vecL(idx_d));
end
[U_tau,dU_tau]=geom_tau(tau,prior_tau,data,vecL);
[U_vecL,dU_vecL]=geom_L(vecL,prior_L,data,tau,~reweight);

% test gradients with finite difference
% test tau
fprintf('\nTesting the difference between true gradient of tau in a random direction and its finite-difference estimate...\n');
v_tau=randn(sqdim,1); 
for i=1:6
    h=10^(-i);
    u_tau_fd=(geom_tau(tau+h.*v_tau,prior_tau,data,vecL)-geom_tau(tau-h.*v_tau,prior_tau,data,vecL))./(2*h);
    gv_tau=v_tau'*dU_tau;
    fprintf('The difference between truth and the estimate is %.8f at stepsize %e.\n',abs(u_tau_fd-gv_tau),h);
end

% test vecL
fprintf('\nTesting the difference between true gradient of vecL in a random direction and its finite-difference estimate...\n');
v_vecL=randn(dim,1);
for i=1:6
    h=10^(-i);
    u_vecL_fd=(geom_L(vecL+h.*v_vecL,prior_L,data,tau,~reweight)-geom_L(vecL-h.*v_vecL,prior_L,data,tau,~reweight))./(2*h);
    gv_vecL=v_vecL'*dU_vecL;
    fprintf('The difference between truth and the estimate %.8f at stepsize %e.\n',abs(u_vecL_fd-gv_vecL),h);
end

