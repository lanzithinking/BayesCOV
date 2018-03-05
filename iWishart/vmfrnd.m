% generate sample(s) from von Mises-Fisher distribution

function vecL = vmfrnd(kappa,mu,n)
if ~exist('n','var')
    n=1;
end
if exist('tmp','dir')~=7
    mkdir('tmp');
end
% set parameters
csvwrite('./tmp/mu.csv',mu);
csvwrite('./tmp/setting.csv',[n,kappa]);
% use R function to get random samples
system('/usr/local/bin/R CMD BATCH rvmf_L.R ./tmp/output4debug.txt');
% retrieve samples
vecL=csvread('./tmp/rv_rmf.csv');

end
