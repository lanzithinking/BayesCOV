% Energy function for cholesky factor L %

function [U,dU]=geom_L(vecL,prior,data,tau,adj_metric)
U=[];dU=[];
if ~exist('adj_metric','var')
    adj_metric=false;
end

d=length(tau);
isigma=exp(-tau);
y_std=(data.y'-data.mu0).*isigma;

L=ivech(vecL,'row');
if contains(prior.choice,'bartlett')
    L=L(d:-1:1,d:-1:1); % becomes upper triangular matrix
end
diagL=diag(L);
invLy_std=L\y_std;
loglik=-data.size*sum(log(abs(diagL)))-sum(invLy_std(:).^2)./2;

switch lower(prior.choice)
    case 'bartlett'
        isigmaitU=isigma.*inv(L)'; % should broadcast in new Matlab
        PsisigmaitU=prior.Psi*isigmaitU;
        dg0=sum(PsisigmaitU.*isigmaitU,2);
        logpri=((1:d)-(prior.nu+d+1))*log(abs(diagL))-0.5*sum(dg0);
    case 'sqdir'
        logpri=(2*prior.alpha-1)'*log(abs(vecL));
    case 'vmf'
        logpri=prior.kappa*sum(diagL);
    case 'bing'
        logpri=prior.zeta*sum(diagL.^2);
    otherwise
        error('Prior not defined!');
end

if adj_metric
    logdetG=-sum(log(abs(diagL)));
else
    logdetG=0;
end

U=-(loglik+logpri+logdetG);
if nargout==2
    dl=-diag(data.size./diagL)+L'\invLy_std*invLy_std';
    if contains(prior.choice,'bartlett')
        dl=dl(d:-1:1,d:-1:1);
    end
    dloglik=vech(tril(dl),'row');
    switch lower(prior.choice)
        case 'bartlett'
            dg=diag(((1:d)'-(prior.nu+d+1))./diagL)+triu(L'\PsisigmaitU'*isigmaitU);
            dlogpri=vech(dg(d:-1:1,d:-1:1),'row');
        case 'sqdir'
            dlogpri=(2*prior.alpha-1)./vecL;
        case 'vmf'
            dlogpri=prior.kappa.*vech(eye(d),'row');
        case 'bing'
            dlogpri=2*prior.zeta.*vech(diag(diagL),'row');
        otherwise
            error('Prior not defined!');
    end
    dlogdetG=zeros(length(vecL),1);
    if adj_metric
        dlogdetG((1:d).*((1:d)+1)./2)=-1./vecL((1:d).*((1:d)+1)./2); % row order
    end
    dU=-(dloglik+dlogpri+dlogdetG);
end

end