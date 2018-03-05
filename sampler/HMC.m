%%%% This is generic HMC algorithm %%%%

function [q,u,du,acpt] = HMC(q_cur,u_cur,du_cur,U,eps,L)

q = q_cur;

% sample velocity
v = randn(length(q),1);

% current energy
E_cur = u_cur + (v'*v)/2;

randL=ceil(rand*L);

% forward half step of velocity
v = v - eps/2.*du_cur;
for i=1:randL
    
    % full step evolution of position
    q = q + eps.*v;
    
    [u,du] = U(q);
    % backward half step of velocity
    if(i~=randL)
        v = v - eps.*du;
    end
    
end
v = v - eps/2.*du;

% new energy
E_prp = u + (v'*v)/2;

% Accept according to ratio
logRatio = -E_prp + E_cur;

if (isfinite(logRatio) && (log(rand) < min([0,logRatio])))
    acpt = 1;
else
    q = q_cur; u = u_cur; du = du_cur;
    acpt = 0;
end


end