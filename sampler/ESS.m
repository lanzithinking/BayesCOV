%%%% Elliptical Slice Sampler algorithm by Murray et~al (2010) %%%%

% inputs:
%   q0: initial state of the parameters
%   l0: initial log-likelihood
%   rnd_pri: random sample from Gaussian prior N(0,C)
%   loglik: log-likelihood function of q
% outputs:
%   q: new state of the parameter following N(q;0,C)*lik
%   l: new log-likelihood


function [q,l] = ESS(q0,l0,rnd_pri,loglik)

% initialization
% q = q0; l = l0;

% sample velocity
v = rnd_pri(); % v ~ N(0,C)
% (q,v) now defines a slice

% log-likelihood threshhold (defines a slice)
logy = l0 + log(rand);

% draw a initial proposal, also defines a bracket
t = 2*pi*rand;
t_min = t-2*pi; t_max = t;

% repeat slice procedure until a proposal is accepted (land on the slice)
while 1
    q = q0.*cos(t) + v.*sin(t);
    l = loglik(q);
%     fprintf('angle between q and v: %.2f\n',acos(q'*v));
    if l>logy
        break;
    end
    % shrink the bracket and try a new point
    if t<0
        t_min = t;
    else
        t_max = t;
    end
    t = t_min+(t_max-t_min)*rand;
%     fprintf('bracket size: %.2f\n',abs(t_min-t_max));
end


end