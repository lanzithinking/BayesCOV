% This is to generate simulated data of a periodic process
% for multiple trials

function [t,y]=generate_data(M,N,D,seedNO)
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

file_name=['periodic_multi_M',num2str(M),'_N',num2str(N),'_D',num2str(D),'.mat'];
if exist(['./',file_name],'file')
    % load data
    load(['./',file_name],'t','y');
else
    % or generate data
    [t,y,~,~,~]=periodic_multiproc(M,N,D,seedNO,true,true);
end

t=t';

end