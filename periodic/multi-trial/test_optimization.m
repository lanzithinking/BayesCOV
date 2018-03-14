% test optimization of initial points for MCMC
% for multiple-trials

clear;
addpath('../../util/','../../optimizer/');
% Random Numbers...
seedNO=2017;
seed=RandStream('mt19937ar','Seed',seedNO);
RandStream.setGlobalStream(seed);

% setting for simulation
M=100; N=200; sqdim=2;
% load or simulate data
[t,y]=generate_data(M,N,sqdim,seedNO);
dim=sqdim*(sqdim+1)/2;
% constant matrix
comm_ND=sparse(vecperm(N,sqdim)');

% mean in uvGP for vecL
M_vecL=ones(N,1)*vech(eye(sqdim),'row')'; % (N,dim)
% % M_=chol(corr(reshape(y,[],sqdim)),'lower');
% M_=corr(reshape(y,[],sqdim));
% M_=sign(M_).*double(abs(M_)>0.9);M_(1:sqdim+1:end)=1;
% M_vecL=ones(N,1)*vech(M_,'row')'; % (N,dim)
% parameters of kernel
s=2; % smoothness
dist_t=pdist2(t,t,@(XI,XJ)abs(bsxfun(@minus,XI,XJ)).^s);
jit=1e-5.*speye(N);
ker.s=s;ker.dist_t=dist_t;ker.jit=jit;
% specify (hyper)-priors
a=ones(1,3); b=[1,1e-3,2e-1]; % (a,b) in inv-gamma priors for gamma_*, * = mu, tau(log-sigma), L
m=zeros(1,3); V=[1,5e-1,2e-1]; % (m,V) in (log) normal priors for eta_*, (eta=log-rho), * = mu, tau(log-sigma), L

% optimize (gamma, eta), mu, tau, L for initializing MCMC.

% initializatioin
gamma=1./gamrnd(a,1./b);
eta=normrnd(m,sqrt(V));
mu=mvnrnd(zeros(1,N),gamma(1).*exp(-.5.*dist_t.*exp(-s.*eta(1)))+jit,sqdim)'; % (N,sqdim)
tau=mvnrnd(zeros(1,N),gamma(2).*exp(-.5.*dist_t.*exp(-s.*eta(2)))+jit,sqdim)'; % (N,sqdim)
vecL=mvnrnd(M_vecL',gamma(3).*exp(-.5.*dist_t.*exp(-s.*eta(3)))+jit)'; % (N,dim)
for d=1:sqdim
    idx_d=1+d*(d-1)/2:d*(d+1)/2;
    vecL(:,idx_d)=vecL(:,idx_d)./sqrt(sum(vecL(:,idx_d).^2,2));
end
vecL(:,1)=1;

% optimize initial location
[gamma,eta,mu,tau,vecL,objf]=opt4ini(gamma,eta,mu,tau,vecL,y,comm_ND,M_vecL,ker,a,b,m,V,[],20);
