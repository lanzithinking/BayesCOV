% this is to test the calculation of gradients

clear;
addpath('../util/');
% Random Numbers...
seedNO=2017;
seed = RandStream('mt19937ar','Seed',seedNO);
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
% parameters of kernel
s=2; % smoothness
dist_t=pdist2(t,t,@(XI,XJ)abs(bsxfun(@minus,XI,XJ)).^s);
jit=1e-7.*speye(N);
ker.s=s;ker.dist_t=dist_t;ker.jit=jit;
% specify (hyper)-priors
a=ones(1,3); b=[1,0.1^2,0.1]; % (a,b) in inv-gamma priors for gamma_*, * = mu, tau(log-sigma), L
m=zeros(1,3); V=ones(1,3); % (m,V) in (log) normal priors for eta_*, (eta=log-rho), * = mu, tau(log-sigma), L

% calculation setting
reweight=false;

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


Rho_mat=zeros(sqdim,sqdim,N);
for i=1:sqdim
    for j=1:sqdim
        ind_i=1+i*(i-1)/2:i*(i+1)/2;
        ind_j=1+j*(j-1)/2:j*(j+1)/2;
        min_ij=min([i,j]);
        Rho_mat(i,j,:)=sum(vecL(:,ind_i(1:min_ij)).*vecL(:,ind_j(1:min_ij)),2);
    end
end
[U_tau,dU_tau]=geom_tau(tau,mu,y,Rho_mat,comm_ND,ker,gamma(2),eta(2));
[U_vecL,dU_vecL]=geom_L(vecL,mu,tau,y,comm_ND,M_vecL,ker,gamma(3),eta(3),[0,1],~reweight);

% test gradients with finite difference
% test tau
fprintf('\nTesting the difference between true gradient of tau in a random direction and its finite-difference estimate...\n');
v_tau=randn(N,sqdim);
for i=1:6
    h=10^(-i);
    u_tau_fd=(geom_tau(tau+h.*v_tau,mu,y,Rho_mat,comm_ND,ker,gamma(2),eta(2))...
              -geom_tau(tau-h.*v_tau,mu,y,Rho_mat,comm_ND,ker,gamma(2),eta(2)))./(2*h);
    gv_tau=dU_tau(:)'*v_tau(:);
    fprintf('The difference between truth and the estimate %.8f at stepsize %e.\n',abs(u_tau_fd-gv_tau),h);
end

% test vecL
fprintf('\nTesting the difference between true gradient of vecL in a random direction and its finite-difference estimate...\n');
v_vecL=randn(N,dim);
for i=1:6
    h=10^(-i);
    u_vecL_fd=(geom_L(vecL+h.*v_vecL,mu,tau,y,comm_ND,M_vecL,ker,gamma(3),eta(3),0,~reweight)...
              -geom_L(vecL-h.*v_vecL,mu,tau,y,comm_ND,M_vecL,ker,gamma(3),eta(3),0,~reweight))./(2*h);
    gv_vecL=dU_vecL(:)'*v_vecL(:);
    fprintf('The difference between truth and the estimate %.8f at stepsize %e.\n',abs(u_vecL_fd-gv_vecL),h);
end