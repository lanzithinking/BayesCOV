% generate sample(s) from Bingham distribution

function vecL = bingrnd(zeta,mu,n)
if ~exist('n','var')
    n=1;
end
if exist('tmp','dir')~=7
    mkdir('tmp');
end
% set parameters
csvwrite('./tmp/mu.csv',mu);
csvwrite('./tmp/setting.csv',[n,zeta]);
% use R function to get random samples
system('/usr/local/bin/R CMD BATCH rbing_L.R ./tmp/output4debug.txt');
% retrieve samples
vecL=csvread('./tmp/rv_bing.csv');

end