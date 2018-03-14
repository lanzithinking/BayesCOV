% This is to generate simulated data of a periodic process

function [t,y]=generate_data(seedNO)
if ~exist('seed','var')
    seedNO=2017;
end

if exist('./periodic_data.mat','file')
    % load data
    load('./periodic_data.mat','t','y');
else
    % or generate data
    % Random Numbers...
    seed = RandStream('mt19937ar','Seed',seedNO);
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
        for j=1:p
            L(i,j,:)=(-1)^(i+j).*sin(i.*t).*cos(j.*t);
            S(i,j)=abs(i-j)+1;
        end
    end
    S=S+tril(S,-1)';
    % form Sigma
    Sigma=squeeze(sum(reshape(L,[p,1,p,N]).*reshape(L,[1,p,p,N]),3))./S; % (p,p,N)

    % generate data
    y=mvnrnd(mu',Sigma)'; % (p,N)

    % save data to file
    save('./periodic_data.mat');
end

end