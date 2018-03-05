%%%% Robust Spherical HMC in Cartesian coordinate %%%%
% for distributions defined on sphere %

% inputs:
%   q_cur: initial state of q with |q|=1
%   u_cur, du_cur: initial potential energy and its gradient
%   U: =-log(density(q)), potential function of q, or its gradient
%   h: step size
%   L: number of leapfrogs
% outputs:
%   q: new state
%   u, du: new potential energy and its gradient
%   acpt: proposal acceptance indicator

function [q, u, du, acpt] = RSphHMC(q_cur, u_cur, du_cur, U, h, L)
if nargin<6
    L=1;
end

% initialization
q = q_cur; D = length(q);
u = u_cur; du = du_cur;

% sample velocity
v = randn(D,1);
v = v - q.*(q'*v); % project to tangent plane

% natural gradient
g = du - q.*(q'*du);

% accumulate the power of force
pow = h/2*(v'*du);
% pow = h/2*(v(1:end-1)'*du(1:end-1));

% current energy
E_cur = u - h^2/8*(du'*g);

% set random integration steps
% randL=ceil(rand*L);
randL=L;

% alternate updating for position and velocity
for l = 1:randL
    % forward half step of velocity
    v = v - h/2.*g;
    
    % full step evolution on sphere along great circle
    q0 = q; v_nom = sqrt(v'*v);
    cosvt = cos(h*v_nom); sinvt = sin(h*v_nom);
    q = q0.*cosvt + v/v_nom.*sinvt;
    v = -q0*v_nom.*sinvt + v.*cosvt;
%     v_nom = sqrt(v'*v);
%     gcirc = (q+1i*v./v_nom).*exp(-1i*v_nom*h);
%     q = real(gcirc); v = imag(gcirc).*v_nom;
    
    % update geometry
    [u,du] = U(q);
%     du = U(q,1);
    g = du - q.*(q'*du);
    
    % backward full step of velocity
    v = v - h/2.*g;
    
    % accumulate the power of force
    if l~=randL
        pow = pow + h*(du'*v);
%         pow = pow + h*(v(1:end-1)'*du(1:end-1));
    end
    
%     % calibrate direction possibly deviated by error accumulation
%     if abs(q'*v)>1e-10
%         v = v - q*(q'*v);
%         disp('Direction calibrated!');
%     end
    
    if q'*q_cur<0
%         disp(['It takes ', num2str(l),' steps to cross 2 orthants!']);
        break;
    end
    
end

% accumulate the power of force
pow = pow + h/2*(du'*v);
% pow = pow + h/2*(v(1:end-1)'*du(1:end-1));

% new energy
% u = U(q);
E_prp = u - h^2/8*(du'*g);

% log of Metropolis ratio
logRatio = -E_prp + E_cur + pow;

% accept or reject the proposal at end of trajectory
if isfinite(logRatio) && (log(rand) < min([0,logRatio]))
    acpt = 1;
else
    q = q_cur; u = u_cur; du = du_cur;
    acpt = 0;
end

end