%%%% Spherical HMC on spheres %%%%

% inputs:
%   q_cur: initial state (matrix) of q, each row satisfying |q_i|=1
%   u_cur, du_cur: initial potential energy and its gradient
%   U: =-log(density(q)), potential function of q, or its gradient
%   h: step size
%   L: number of leapfrogs
% outputs:
%   q: new state
%   u, du: new potential energy and its gradient
%   acpt: proposal acceptance indicator

function [q, u, du, acpt] = SphHMCs(q_cur, u_cur, du_cur, U, h, L)
if nargin<6
    L=1;
end

% initialization
q = q_cur; D = size(q);
u = u_cur; du = du_cur;

% sample velocity
v = randn(D); % same size of state
v = v - bsxfun(@times,sum(v.*q,2),q); % project to tangent plane

% current energy
E_cur = u + sum(sum(v.^2))/2;

% set random integration steps
randL=ceil(rand*L);

% forward half step of velocity
g = du - bsxfun(@times,sum(du.*q,2),q);
v = v - h/2.*g;
for l = 1:randL
    
    % full step evolution on sphere along great circle
    q0 = q; v_nom = sqrt(sum(v.^2,2));
    cosvt = cos(h.*v_nom); sinvt = sin(h.*v_nom);
    q = bsxfun(@times,q0,cosvt) + bsxfun(@times,v,sinvt./v_nom);
    v = -bsxfun(@times,q0,v_nom.*sinvt) + bsxfun(@times,v,cosvt);
    
%     rot = bsxfun(@times,q+1i*bsxfun(@rdivide,v,v_nom),exp(1i*(v_nom.*h)));
%     q = real(rot); v = bsxfun(@times,imag(rot),v_nom);
    
    % backward full step of velocity
    du = U(q,1); g = du - bsxfun(@times,sum(du.*q,2),q);
    if l~=randL
        v = v - h.*g;
    end
    
    % calibrate direction possibly deviated by error accumulation
    if any(abs(sum(q.*v,2))>1e-6)
        v = v - bsxfun(@times,sum(v.*q,2),q);
        disp('Direction calibrated!');
    end
    
end
% backward last half step of velocity
v = v - h/2.*g;

% new energy
u = U(q);
E_prp = u + sum(sum(v.^2))/2;

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