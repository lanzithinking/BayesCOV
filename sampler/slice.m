%%%% Slice sampling (1d) algorithm by Neal (2003) %%%%

% inputs:
%   q0: initial state of the parameters
%   l0: initial log-density
%   logf: log-density function of q
%   w: estimate size of a slice
%   m: integer limiting th size of a slice to mw
% outputs:
%   q: new state of the parameter following N(q;0,C)*lik
%   l: new log-density


function [q,l] = slice(q0,l0,logf,w,m,bdy)
if ~exist('w','var')
    w=1;
end
if ~exist('m','var')
    m=5;
end
if ~exist('bdy','var')
    bdy=[-Inf,+Inf];
end

% log-denisty threshhold (defines a slice)
logy = l0 + log(rand);

% stepp out to obtain the [L,R] range
L = q0 - w.*rand;
R = L + w;
J = floor(m.*rand);
K = (m-1) - J;

% make sure [L,R] is within boundary
L = max([L,bdy(1)]);
R = min([R,bdy(2)]);

while 1
    if J<=0 || logy>=logf(L)
        break
    end
    L = L - w;
    L = max([L,bdy(1)]);
    J = J - 1;
end

while 1
    if K<=0 || logy>=logf(R)
        break;
    end
    R = R + w;
    R = min([R,bdy(2)]);
    K = K - 1;
end


% skrink to obtain a sample
while 1
    q = L + rand.*(R-L);
    l = logf(q);
    if l>logy
        break;
    end
    % shrink the bracket and try a new point
    if q<q0
        L = q;
    else
        R = q;
    end
end

end