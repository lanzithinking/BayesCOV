% Energy function for log-variances tau %

function [U,dU]=geom_tau(tau,prior,data,vecL)
U=[];dU=[];

d=length(tau);
isigma=exp(-tau);
y_std=(data.y'-data.mu0).*isigma;

L=ivech(vecL,'row');
if contains(prior.choice,'bartlett')
    L=L(d:-1:1,d:-1:1); % becomes upper triangular matrix
end
invLy_std=L\y_std;
loglik=-data.size*sum(tau)-sum(invLy_std(:).^2)./2;

switch lower(prior.choice)
    case 'bartlett'
        isigmaitU=isigma.*inv(L)'; % should broadcast in new Matlab
        PsisigmaitU=prior.Psi*isigmaitU;
        dg=sum(PsisigmaitU.*isigmaitU,2);
        logpri=((1:d)-(prior.nu+d))*tau-0.5*sum(dg);
    case 'normal'
        invCtau=(tau-prior.mean)./prior.var;
        logpri=-(tau-prior.mean)'*invCtau./2;
    otherwise
        error('Prior not defined!');
end

U=-(loglik+logpri);
if nargout==2
    dloglik=-data.size+sum(y_std.*(L'\invLy_std),2);
    switch lower(prior.choice)
        case 'bartlett'
            dlogpri=(1:d)'-(prior.nu+d)+dg;
        case 'normal'
            dlogpri=-invCtau;
        otherwise
            error('Prior not defined!');
    end
    dU=-(dloglik+dlogpri);
end

end