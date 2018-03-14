% test optimization of initial points for MCMC

clear;
addpath('../util/','../optimizer/');
% Random Numbers...
seedNO=2017;
seed=RandStream('mt19937ar','Seed',seedNO);
RandStream.setGlobalStream(seed);

% load or simulate data
[t,y]=generate_data(seedNO);
N=length(t);
sqdim=size(y,1);
dim=sqdim*(sqdim+1)/2;
t=t'; y=y';
% constant matrix
comm_ND=sparse(vecperm(N,sqdim)');

% mean in uvGP for vecL
M_vecL=ones(N,1)*vech(eye(sqdim),'row')'; % (N,dim)
% M=chol(corr(y),'lower');
% M=eye(sqdim); M(2,1)=.1;
% M_vecL=ones(N,1)*vech(M,'row')'; % (N,dim)
% parameters of kernel
s=2; % smoothness
dist_t=pdist2(t,t,@(XI,XJ)abs(bsxfun(@minus,XI,XJ)).^s);
jit=1e-5.*speye(N);
ker.s=s;ker.dist_t=dist_t;ker.jit=jit;
% specify (hyper)-priors
a=ones(1,3); b=[1,1e-2,2e-1]; % (a,b) in inv-gamma priors for gamma_*, * = mu, tau(log-sigma), L
m=[0,0,0]; V=[1,1,1e-2]; % (m,V) in (log) normal priors for eta_*, (eta=log-rho), * = mu, tau(log-sigma), L

% need to sample mu, tau, L, gamma, eta
% using Gibbs, ESS, Dlt-SphHMC, Gibbs, slice sampling methods resp.


% initializatioin
gamma=1./gamrnd(a,1./b);
eta=normrnd(m,sqrt(V));
mu=mvnrnd(zeros(1,N),gamma(1).*(exp(-.5.*dist_t.*exp(-s.*eta(1)))+jit),sqdim)'; % (N,sqdim)
tau=mvnrnd(zeros(1,N),gamma(2).*(exp(-.5.*dist_t.*exp(-s.*eta(2)))+jit),sqdim)'; % (N,sqdim)
vecL=mvnrnd(M_vecL',gamma(3).*(exp(-.5.*dist_t.*exp(-s.*eta(3)))+jit))'; % (N,dim)
for d=1:sqdim
    idx_d=1+d*(d-1)/2:d*(d+1)/2;
    vecL(:,idx_d)=vecL(:,idx_d)./sqrt(sum(vecL(:,idx_d).^2,2));
end
vecL(:,1)=1;

% optimize initial location
[gamma,eta,mu,tau,vecL,objf]=opt4ini(gamma,eta,mu,tau,vecL,y,comm_ND,M_vecL,ker,a,b,m,V,[],20);
