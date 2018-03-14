% This is to generate simulated data of multiple periodic process

function [t,y,mu,Sigma,Rho]=periodic_multiproc(M,N,D,seedNO,SAVE,PLOT)
if ~exist('M','var')
    M=10; % number of trials
end
if ~exist('N','var')
    N=100; % discretization size
end
if ~exist('D','var')
    D=2; % data dimension
end
if ~exist('seed','var')
    seedNO=2017;
end
if ~exist('SAVE','var')
    SAVE=false;
end
if ~exist('PLOT','var')
    PLOT=false;
end

% addpath
addpath('../../util/');
% Random Numbers...
seed = RandStream('mt19937ar','Seed',seedNO);
RandStream.setGlobalStream(seed);

% parameters setting
% t=1/N:1/N:1;
t=linspace(0,2,N+1); t=t(2:end);

% set mu
mu=zeros(D,N);
for i=1:D
    mu(i,:)=sin(i.*t.*pi/D);
end

% set L
L=zeros(D,D,N); S=zeros(D);
for i=1:D
    for j=1:i
        L(i,j,:)=(-1)^(i+j).*sin(i.*t.*pi/D).*cos(j.*t.*pi/D);%+i/D*j/D;
        S(i,j)=abs(i-j)+1;
    end
end
S=S+tril(S,-1)';
% form Sigma
Sigma=squeeze(sum(reshape(L,[D,1,D,N]).*reshape(L,[1,D,D,N]),3))./S; % (D,D,N)

% get Rho
Rho=zeros(D,D,N);
for n=1:N
    Rho(:,:,n)=diag(sqrt(diag(Sigma(:,:,n))))\Sigma(:,:,n)/diag(sqrt(diag(Sigma(:,:,n))));
end

% % set Sigma to unit-variance to focus on correlation
% Sigma=Rho;

% generate data
y=zeros(M,N,D);
for m=1:M
    y(m,:,:)=mvnrnd(mu',Sigma); % (N,D)
end

% save data to file
if SAVE
    save(['./periodic_multi_M',num2str(M),'_N',num2str(N),'_D',num2str(D),'.mat']);
end

% plot data
if PLOT
    addpath('~/Documents/MATLAB/tight_subplot/');
    addpath('~/Documents/MATLAB/columnlegend/');
    % setting
    fig=figure(1); clf(fig);
    set(fig,'pos',[0 800 1000 350]);
    ha=tight_subplot(1,3,[0.1,.06],[.13,.08],[.05,.04]);
    % plot y_1 and mean
    subplot(ha(1));
    plot(t',squeeze(y(1,:,:)),'*','markersize',4.5); hold on;
    set(gca,'colororderindex',1);
    plot(t',mu','linewidth',2);
    xlim([min(t),max(t)]); ylim([-3,3]);
    set(gca,'fontsize',14);
    xlabel('t','fontsize',16); ylabel('y','rot',0,'fontsize',16);
    lgd_str=cell(2*D,1);
    lgd_str(1:D)=cellstr([repmat('y_',D,1),num2str((1:D)')]);
    lgd_str(1+D:end)=cellstr([repmat('\mu_',D,1),num2str((1:D)')]);
    lgd=columnlegend(D,lgd_str,'location','southwest');
    set(lgd,'fontsize',14);
    % plot variance
    subplot(ha(2));
    I=1:D; J=1:D;
    for i=1:length(I)
        plot(t',squeeze(Sigma(I(i),J(i),:)),'displayname',['\sigma^2_{',num2str(I(i)),'}']); hold on;
    end
    xlim([min(t),max(t)]); ylim([-0.2,1]);
    set(gca,'fontsize',14);
    xlabel('t','fontsize',16); ylabel('Variance','fontsize',16);
    lgd=legend('location','southeast');
    set(lgd,'orientation','horizontal','fontsize',14,'box','off');
    % plot correlation
    subplot(ha(3));
    [I,J]=ind_vech2sub(D,setdiff(1:D*(D+1)/2,(1:D).*((1:D)+1)./2),'row');
    for i=1:length(I)
        plot(t',squeeze(Rho(I(i),J(i),:)),'displayname',['\rho_{',num2str(J(i)),',',num2str(I(i)),'}']); hold on;
    end
    xlim([min(t),max(t)]); ylim([-0.6,0.6]);
    set(gca,'fontsize',14);
    xlabel('t','fontsize',16); ylabel('Correlation','fontsize',16);
    lgd=legend('location','southeast');
    set(lgd,'orientation','horizontal','fontsize',14,'box','off');
    % save plot
    if SAVE
        fig.PaperPositionMode = 'auto';
        print(fig,['./periodic_multi_M',num2str(M),'_N',num2str(N),'_D',num2str(D)],'-dpng','-r0');
    end
end

end