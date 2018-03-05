%%%% Spherical HMC in Cartesian coordinate %%%%
% for distributions defined on sphere %

% inputs:
%   q_cur: initial state of q, with |q|=1
%   u_cur, du_cur: initial potential energy and its gradient
%   U: =-log(density(q)), potential function of q, or its gradient
%   h: step size
%   L: number of leapfrogs
% outputs:
%   q: new state
%   u, du: new potential energy and its gradient
%   acpt: proposal acceptance indicator

function [q, u, du, acpt] = cSphHMC(q_cur, u_cur, du_cur, U, h, L)
if nargin<6
    L=1;
end

% initialization
q = q_cur; D = length(q);
u = u_cur; du = du_cur;

% sample velocity
v = randn(D,1);
v = v - q.*(q'*v); % project to tangent plane

% current energy
E_cur = u + (v'*v)/2;

% set random integration steps
randL=ceil(rand*L);

% forward half step of velocity
g = du - q.*(q'*du);
v = v - h/2.*g;
for l = 1:randL
    
    % full step evolution on sphere along great circle
    q0 = q; v_nom = sqrt(v'*v);
    cosvt = cos(h*v_nom); sinvt = sin(h*v_nom);
    q = q0.*cosvt + v/v_nom.*sinvt;
    v = -q0*v_nom.*sinvt + v.*cosvt;
    
    % backward full step of velocity
    du = U(q,1);
    g = du - q.*(q'*du);
    if l~=randL
        v = v - h.*g;
    end
    
%     % calibrate direction possibly deviated by error accumulation
%     if abs(q'*v)>1e-6
%         v = v - q*(q'*v);
%         disp('Direction calibrated!');
%     end
    
end
% backward last half step of velocity
v = v - h/2.*g;

% new energy
u = U(q);
E_prp = u + (v'*v)/2;

% log of Metropolis ratio
logAP = -E_prp + E_cur;

% accept or reject the proposal at end of trajectory
if isfinite(logAP) && (log(rand) < min([0,logAP]))
    acpt = 1;
else
    q = q_cur; u = u_cur; du = du_cur;
    acpt = 0;
end

end