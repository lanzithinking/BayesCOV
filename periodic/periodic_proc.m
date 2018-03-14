% This is to generate simulated data of a periodic process

clear;
addpath('../util/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',2017);
RandStream.setGlobalStream(seed);

% parameters setting
N=200; % discretization size
p=2; % data dimension
% t=1:N;
% t=1/N:1/N:1;
t=linspace(0,2*pi/p,N+1); t=t(2:end);

% set mu
mu=zeros(p,N);
for i=1:p
    mu(i,:)=sin(i.*t);
end

% set L
L=zeros(p,p,N); S=zeros(p);
for i=1:p
    for j=1:i
        L(i,j,:)=(-1)^(i+j).*sin(i.*t).*cos(j.*t);
        S(i,j)=abs(i-j)+1;
    end
end
S=S+tril(S,-1)';
% form Sigma
Sigma=squeeze(sum(reshape(L,[p,1,p,N]).*reshape(L,[1,p,p,N]),3))./S; % (p,p,N)

% get Rho
Rho=zeros(p,p,N);
for n=1:N
    Rho(:,:,n)=diag(sqrt(diag(Sigma(:,:,n))))\Sigma(:,:,n)/diag(sqrt(diag(Sigma(:,:,n))));
end

% generate data
y=mvnrnd(mu',Sigma)'; % (p,N)

% save data to file
save('./periodic_data.mat');

% plot data
addpath('~/Documents/MATLAB/tight_subplot/');
addpath('~/Documents/MATLAB/columnlegend/');
% setting
fig=figure(1); clf(fig);
set(fig,'pos',[0 800 1000 350]);
ha=tight_subplot(1,3,[0.1,.06],[.13,.08],[.05,.04]);
% plot y and mean
subplot(ha(1));
plot(t,y,'*','markersize',4.5); hold on;
set(gca,'colororderindex',1);
plot(t',mu','linewidth',2);
xlim([0,2*pi/p]);
set(gca,'fontsize',14);
xlabel('t','fontsize',16); ylabel('y','rot',0,'fontsize',16);
lgd_str=cell(2*p,1);
lgd_str(1:p)=cellstr([repmat('y_',p,1),num2str((1:p)')]);
lgd_str(1+p:end)=cellstr([repmat('\mu_',p,1),num2str((1:p)')]);
lgd=columnlegend(2,lgd_str,'location','southwest');
set(lgd,'fontsize',14);
% % plot mean
% subplot(ha(2));
% plot(t',mu');
% xlim([0,2*pi/p]);
% set(gca,'fontsize',14);
% xlabel('t','fontsize',16); ylabel('Mean','fontsize',16);
% lgd=legend([repmat('\mu_',p,1),num2str((1:p)')]);
% set(lgd,'location','southwest','fontsize',14);
% plot variance
subplot(ha(2));
% plot(t',reshape(Sigma,[],N)');
% [I,J]=ind_vech2sub(p,1:p*(p+1)/2,'row');
[I,J]=ind_vech2sub(p,[1,3,2],'row');
for i=1:length(I)
    plot(t',squeeze(Sigma(I(i),J(i),:)),'displayname',['\sigma_{',num2str(J(i)),',',num2str(I(i)),'}']); hold on;
end
xlim([0,2*pi/p]); ylim([-0.5,1]);
set(gca,'fontsize',14);
xlabel('t','fontsize',16); ylabel('Variance/Covariance','fontsize',16);
lgd=legend('location','southeast');
set(lgd,'orientation','horizontal','fontsize',14);
% plot correlation
subplot(ha(3));
% plot(t',reshape(Sigma,[],N)');
% [I,J]=ind_vech2sub(p,1:p*(p+1)/2,'row');
[I,J]=ind_vech2sub(p,2,'row');
for i=1:length(I)
    plot(t',squeeze(Rho(I(i),J(i),:)),'displayname',['\rho_{',num2str(J(i)),',',num2str(I(i)),'}']); hold on;
end
xlim([0,2*pi/p]);
set(gca,'fontsize',14);
xlabel('t','fontsize',16); ylabel('Correlation','fontsize',16);
lgd=legend('location','southeast');
set(lgd,'fontsize',14);
% save plot
fig.PaperPositionMode = 'auto';
print(fig,'./periodic_data','-dpng','-r0');