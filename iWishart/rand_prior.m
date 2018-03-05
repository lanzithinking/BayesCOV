% random number generator for prior distributions %

function vecL=rand_prior(prior,tau)
if ~exist('tau','var')
    tau=0;
end

sqdim=(sqrt(8*prior.dim+1)-1)/2;
vecL=zeros(prior.dim,1);

switch lower(prior.choice)
    case 'bartlett'
        isigma=exp(-tau);
        Sigma=iwishrnd(isigma.*prior.Psi.*isigma',prior.nu); % L|sigma
        L=chol(Sigma(end:-1:1,end:-1:1),'lower');
        L=L./sqrt(sum(L.^2,2)); % should broadcast in new Matlab
        vecL=vech(L,'row');
    case 'sqdir'
        for d=1:sqdim
            idx_d=1+d*(d-1)/2:d*(d+1)/2;
            vecL(idx_d)=sqdrchrnd(prior.alpha(idx_d)');
        end
    case 'vmf'
%         for d=1:sqdim
%             idx_d=1+d*(d-1)/2:d*(d+1)/2;
%             mu=zeros(1,d); mu(end)=1;
%             vecL(idx_d)=vmfrnd(prior.kappa,mu,1);
%         end
        vecL=vmfrnd(prior.kappa,mu);
    case 'bing'
%         for d=1:sqdim
%             idx_d=1+d*(d-1)/2:d*(d+1)/2;
%             mu=zeros(1,d); mu(end)=1;
%             vecL(idx_d)=bingrnd(prior.zeta,mu,1);
%         end
        vecL=bingrnd(prior.zeta,mu);
    otherwise
        error('Prior not defined!');
end


end