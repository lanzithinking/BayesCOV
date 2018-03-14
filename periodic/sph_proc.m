% This is to experiment with spherical processes


clear;
addpath('../util/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% parameters setting
N=100; % discretization size
% t=1:N;
t=1/N:1/N:1;
D=3; % vector process dimension
% hyper-parameters
a_theta=1; b_theta=0.1;
sigma_theta=10;

% sample GP parameters
tau2_theta=1./gamrnd(a_theta,1/b_theta);
rho_theta=exp(normrnd(0,sigma_theta));

% form kernel matrix
K=tau2_theta.*exp(-.5*pdist2(t',t','squaredeuclidean')./rho_theta)+1e-7.*eye(N);
% sample Z from Gaussian process
% Z=mvnrnd(zeros(1,N),K,D*(D+1)/2); % (D*(D+1)/2,N)
Z=mvnrnd(vech(eye(D),'row')*ones(1,N),K); % (D*(D+1)/2,N)

% normalize GP for L
L=zeros(D,D,N);
vecL=zeros(D*(D+1)/2,N);
for d=1:D
    idx_d=1+d*(d-1)/2:d*(d+1)/2;
    % generate spherical process
    % normalized GP
    l=Z(idx_d,:)./sqrt(sum(Z(idx_d,:).^2,2));
    vecL(idx_d,:)=l;
    L(d,1:d,:)=l;
end

% form correlation process
Rho=squeeze(sum(reshape(L,[D,1,D,N]).*reshape(L,[1,D,D,N]),3)); %(D,D,N)

% plot results
addpath('~/Documents/MATLAB/tight_subplot/');
% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 350]);
ha=tight_subplot(1,3,[0.1,.07],[.13,.07],[.06,.04]);
% plot Z
subplot(ha(1));
plot(t',Z');
set(gca,'fontsize',14);
xlabel('t','fontsize',16); ylabel('Latent GPs','fontsize',16);
title('${\bf Z}_t$','interpreter','latex','fontsize',16);
% plot Rho
subplot(ha(2));
plot(t',vecL');
set(gca,'fontsize',14);
xlabel('t','fontsize',16); ylabel('Normalized GPs','fontsize',16);
title('${\bf L}_t$','interpreter','latex','fontsize',16);
% plot Rho
subplot(ha(3));
plot(t',reshape(Rho,[],N)');
set(gca,'fontsize',14);
xlabel('t','fontsize',16); ylabel('Correlations','fontsize',16);
title('${\bf P}_t$','interpreter','latex','fontsize',16);
% save plot
fig.PaperPositionMode = 'auto';
print(fig,'./sph_proc','-dpng','-r0');
