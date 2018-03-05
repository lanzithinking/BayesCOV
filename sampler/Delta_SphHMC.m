%%%% ---- Spherical Hamiltonian Monte Carlo ---- %%%%
%%%% Shiwei Lan                                  %%%%
%%%% CMS, CalTech                                %%%%
%%%% slan@caltech.edu; lanzithinking@outlook.com %%%%
%%%% Copyright @ 2016                            %%%%
%%%% ---- Spherical Hamiltonian Monte Carlo ---- %%%%

classdef Delta_SphHMC
    % Spherical Hamiltonian Monte Carlo is a MCMC algorithm for sampling
    % distributions defined on sphere. Standard algorithm needs hand tuning
    % the step size and the number of leapfrog steps; while adaptive
    % algorithm introduces automatic tuning so users do not have to tune.
    properties
        r=1; % radius of the sphere
        r_xt; % extended radius of the sphere (repeat r to match size of q)
        q; dim; % current state and its dimension
        sqdim; % the dimension of q in square matrix
        U; dU; % current potential and its gradient
        geom; % the function to calculate geometry (potential and gradient)
        h=1; L=20; % the step size and the number of leapfrogs (standard)
        alg='standard'; % the name of algorithm
        h_adpt; % struct for adaptation of step size
        reweight=true; % flag of reweight
    end
    methods
        % constructor
        function obj=Delta_SphHMC(q,geom,h,L,alg,reweight)
            obj.q=q; obj.dim=length(q);
            obj.sqdim=ceil((sqrt(1+8*obj.dim)-1)/2);
            obj.r=ones(obj.sqdim,1);
            obj.r_xt=ones(obj.dim,1);
            for d=1:obj.sqdim
                idx_d=1+d*(d-1)/2:d*(d+1)/2;
                obj.r(d)=sqrt(q(idx_d)'*q(idx_d));
                obj.r_xt(idx_d)=obj.r(d);
            end
            if exist('reweight','var')
                obj.reweight=reweight;
            end
            obj.geom=@(q)geom(q,~obj.reweight);
            [obj.U,obj.dU]=obj.geom(q);
            obj.h=h; obj.L=L;
            if exist('alg','var')
                obj.alg=alg;
            end
            if contains(obj.alg,{'adp','adpt','adaptive'},'IgnoreCase',true)
                obj.alg='adaptive';
%                 h_adpt.h=obj.init_h();
                h_adpt.h=obj.h;
                h_adpt.mu=log(10*h_adpt.h);
                h_adpt.loghn=0;
                h_adpt.An=0;
                % constants' setting
                h_adpt.gamma=0.05;
                h_adpt.n0=10;
                h_adpt.kappa=0.75;
                h_adpt.a0=0.65;
                obj.h_adpt=h_adpt;
            end
        end
        
        % projection a vector v onto each of (d-1) sphere S^(i-1)
        function v=mult_proj(obj,v,q,r)
            if nargin<3
                q=obj.q; r=obj.r;
            end
            for d=1:obj.sqdim
                idx_d=1+d*(d-1)/2:d*(d+1)/2;
                v(idx_d)=v(idx_d)-q(idx_d).*(q(idx_d)'*v(idx_d))./r(d)^2;
            end
        end
        
        % resample v
        function v=resample_aux(obj,q,r)
            if nargin<2
                q=obj.q; r=obj.r;
            end
            v=obj.mult_proj(randn(obj.dim,1),q,r);
        end
        
        % 1 step of leapfrog
        function [q,v,U,dU]=leapfrog(obj,q,v,dU,h,r)
            if nargin<5
                h=obj.h; r=obj.r;
            elseif nargin<6
                r=obj.r;
            end
            
            % natural gradient
            ng=obj.mult_proj(dU,q,r);
            
            % forward half step of velocity
            v=v -h/2.*ng;
            
            % full step evolution on sphere along great circles
            v_nom = zeros(obj.dim,1);
            for d=2:obj.sqdim
                idx_d=1+d*(d-1)/2:d*(d+1)/2;
                v_nom(idx_d)=sqrt(v(idx_d)'*v(idx_d))./r(d);
            end
            q0 = q;
            cosvt = cos(h.*v_nom); sinvt = sin(h.*v_nom);
            q = q0.*cosvt + v./[1;v_nom(2:end)].*sinvt;
            v = -q0.*v_nom.*sinvt + v.*cosvt;
%             gcirc = (q+1i*v./[1;v_nom(2:end)]).*exp(-1i*v_nom.*h);
%             q = real(gcirc); v = imag(gcirc).*v_nom;
            
            % update geometry
            [U,dU] = obj.geom(q);
            ng=obj.mult_proj(dU,q,r);
            
            % backward half step of velocity
            v=v -h/2.*ng;
            
            % adjust the state when necessary
            if abs(q'*v)>1e-6
                [q,v]=obj.adj2S(q,v,r);
            end
        end
        
        % adjust the position to stay on the sphere
        function [q,v]=adj2S(obj,q,v,r)
            if nargin<4
                r=obj.r;
            end
            for d=1:obj.sqdim
                idx_d=1+d*(d-1)/2:d*(d+1)/2;
                q(idx_d)=q(idx_d)./norm(q(idx_d));
                v(idx_d)=v(idx_d)-q(idx_d).*(q(idx_d)'*v(idx_d));
                q(idx_d)=q(idx_d).*r(d);
            end
        end
        
         % standard multiple (Delta) spherical HMC
        function [obj,acpt_ind]=std_Dlt_SphHMC(obj,rndL,h,L,r)
            if nargin<2
                rndL=false; h=obj.h; L=obj.L; r=obj.r;
            elseif nargin<3
                h=obj.h; L=obj.L; r=obj.r;
            elseif nargin<4
                L=obj.L; r=obj.r;
            elseif nargin<5
                r=obj.r;
            end
            % initialization
            q=obj.q; dU=obj.dU;
            v=obj.resample_aux(q,r);
            
            % current energy
            E_cur=obj.U+(v'*v)./2;
            
            % randomize the number of integration steps if asked
            if rndL
                L=ceil(rand*L);
            end
            
            % leapfrogs
            for l=1:L
                [q,v,U,dU]=obj.leapfrog(q,v,dU,h,r);
            end
            
            % proposed energy
            E_prp=U+(v'*v)./2;
            
            % log of Metropolis ratio
            logAP = -E_prp + E_cur;
            % accept or reject the proposal at end of trajectory
            if isfinite(logAP) && (log(rand) < min([0,logAP]))
                obj.q=q; obj.U=U; obj.dU=dU;
                acpt_ind = true;
            else
                acpt_ind = false;
            end
                
        end
        
        % find reasonable initial step size
        function h=init_h(obj)
            % initialization
            q0=obj.q; U0=obj.U; dU0=obj.dU;
            v0=obj.resample_aux();
            E_cur=U0+(v0'*v0)./2;
            h=1;
            [q,v,U,dU]=obj.leapfrog(q0,v0,dU0,h);
            E_prp=U+(v'*v)./2;
            logAP = -E_prp + E_cur;
            a=2*(exp(logAP)>0.5)-1;
            while a*logAP>-a*log(2)
                h=h*2^a;
                [q,v,U,dU]=obj.leapfrog(q0,v0,dU0,h);
                E_prp=U+(v'*v)./2;
                logAP = -E_prp + E_cur;
            end
        end
        
        % dual-averaging to adapt step size
        function hn_adpt=dual_avg(obj,iter,an)
            hn_adpt=obj.h_adpt;
            hn_adpt.An=(1-1./(iter+hn_adpt.n0))*hn_adpt.An + (hn_adpt.a0-an)./(iter+hn_adpt.n0);
            logh=hn_adpt.mu - sqrt(iter)/hn_adpt.gamma.*hn_adpt.An;
            hn_adpt.loghn=iter^(-hn_adpt.kappa).*logh + (1-iter^(-hn_adpt.kappa)).*hn_adpt.loghn;
            hn_adpt.h=exp(logh);
        end
        
        % adaptive Spherical HMC
        function [obj,acpt_ind]=adp_SphHMC(obj,iter,Nadpt)
            % initialization
            q=obj.q; dU=obj.dU;
            v=obj.resample_aux();
            
            % current energy;
            E_cur=obj.U+(v'*v)./2;
            
            % leapfrogs
            for l=1:obj.L
%             for l=1:ceil(1/obj.h_adpt.h)
                [q,v,U,dU]=obj.leapfrog(q,v,dU,obj.h_adpt.h);
                % two orthants' rule
                if q'*obj.q<0
%                 % symmetric stochastic rule
%                 if rand <= (1-q'*obj.q./obj.r^2)/2
                    fprintf('It takes %d steps to exit!\n',l);
                    break;
                end
            end
            
            % proposed energy
            E_prp=U+(v'*v)./2;
            
            % log of Metropolis ratio
            logAP = -E_prp + E_cur;
            % accept or reject the proposal at end of trajectory
            if isfinite(logAP) && (log(rand) < min([0,logAP]))
                obj.q=q; obj.U=U; obj.dU=dU;
                acpt_ind = true;
            else
                acpt_ind = false;
            end
            
            % adapt step size h if needed
            if iter<=Nadpt
                obj.h_adpt=obj.dual_avg(iter,exp(min([0,logAP])));
                fprintf('New step size: %.2f',obj.h_adpt.h);
                fprintf('\tNew averaged step size: %.6f\n',exp(obj.h_adpt.loghn));
            end
            % stop adaptation and freeze the step size
            if iter==Nadpt
                obj.h_adpt.h=exp(obj.h_adpt.loghn);
                fprintf('Adaptation completed; step size freezed at:  %.6f\n',obj.h_adpt.h);
            end
            
        end
        
        % sample function
        function [acpt,time,f_out]=sample(obj,Nsamp,burnrate,thin,SAVE)
            if nargin<4
                thin=1; SAVE=True;
            elseif nargin<5
                SAVE=True;
            end
            % setting
            Niter=Nsamp*thin;
            NBurnIn=floor(Niter*burnrate);
            Niter=Niter+ NBurnIn;
            % allocation to save
            samp=zeros(Nsamp,obj.dim);
            engy=zeros(Niter,1);
            acpt=0; % overall acceptance
            accp=0; % online acceptance
            % run MCMC to collect samples
            fprintf('Running %s Delta Spherical HMC...\n',obj.alg);
            prog=0.05:0.05:1;
            tic;
            for iter=1:Niter

                % display sampleing progress and online acceptance rate
                if ismember(iter,floor(Niter.*prog))
                    fprintf('%.0f%% iterations completed.\n',100*iter/Niter);
                    fprintf('Online acceptance rate: %.2f%%\n', 100*accp/floor(prog(1)*Niter));
                    accp=0;
                end

                % Use Spherical HMC to get samples
                switch obj.alg
                    case 'standard'
                        [obj,acpt_ind]=obj.std_Dlt_SphHMC();
                    case 'adaptive'
                        [obj,acpt_ind]=obj.rbt_Dlt_SphHMC(iter,NBurnIn);
                end
                accp=accp+acpt_ind;

                % burn-in complete
                if(iter==NBurnIn)
                    fprintf('Burn in completed!\n');
                end

                engy(iter) = obj.U;
                % save samples after burn-in
                if(iter>NBurnIn)
                    if mod(iter-NBurnIn,thin) == 0
                        samp(ceil((iter-NBurnIn)/thin),:) = obj.q;
                    end
                    acpt = acpt + acpt_ind;
                end
            end
            % count time
            time=toc;
            % save results
            acpt=acpt/(Niter-NBurnIn);
            % resample if set to reweight
            if obj.reweight
                wt=sum(log(abs(samp(:,(1:obj.sqdim).*((1:obj.sqdim)+1)./2))),2);
                wt=exp(wt-median(wt));
                resamp_idx=datasample(1:Nsamp,Nsamp,'replace',true,'weights',wt);
            else
                wt=ones(Nsamp,1);
                resamp_idx=1:Nsamp;
            end
            % save
            if SAVE
                time_lbl=regexprep(num2str(fix(clock)),'    ','_');
                f_name=[obj.alg,'_Dlt_SphHMC_dim',num2str(obj.dim),'_',time_lbl];
                if exist('result','dir')~=7
                    mkdir('result');
                end
                h=obj.h; L=obj.L;
                save(['result/',f_name,'.mat'],'Nsamp','NBurnIn','thin','h','L','samp','wt','resamp_idx','engy','acpt','time');
                f_out=f_name;
            else
                f_out=struct('samp',samp,'engy',engy,'wt',wt,'resamp_idx',resamp_idx);
            end
        end
    end
end