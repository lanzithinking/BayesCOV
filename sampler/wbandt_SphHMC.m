%%%% ---- Spherical Hamiltonian Monte Carlo ---- %%%%
%%%% Shiwei Lan                                  %%%%
%%%% CMS, CalTech                                %%%%
%%%% slan@caltech.edu; lanzithinking@outlook.com %%%%
%%%% Copyright @ 2017                            %%%%
%%%% ---- Spherical Hamiltonian Monte Carlo ---- %%%%

classdef wbandt_SphHMC
    % Spherical Hamiltonian Monte Carlo is a MCMC algorithm for sampling
    % distributions defined on sphere. Standard algorithm needs hand tuning
    % the step size and the number of leapfrog steps; while adaptive
    % algorithm introduces automatic tuning so users do not have to tune.
    % -- modified for applications of processes (time series) --
    % -- speed up by exploring k-band structure of covariance --
    properties
        r=1; % radius of the sphere
        q; N; D; w; % current state and its series length, dimension and diagonal bandwidth
        U; dU; % current potential and its gradient
        geom; % the function to calculate geometry (potential and gradient)
        h=1; L=20; % the step size and the number of leapfrogs (standard)
        alg='standard'; % the name of algorithm
        h_adpt; % struct for adaptation of step size
        reweight=true; % flag of reweight
    end
    methods
        % constructor
        function obj=wbandt_SphHMC(q,geom,h,L,alg,reweight)
            obj.q=q; [obj.N,obj.D,obj.w]=size(q);
            obj.r=ones(1,obj.D);
            if exist('reweight','var')
                obj.reweight=reweight;
            end
            obj.geom=@(q,geord)geom(q,geord,~obj.reweight);
            [obj.U,obj.dU]=obj.geom(q,[0,1]);
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
        
        % projection a vector v onto each of (d-1) sphere S^w
        function v=mult_proj(obj,v,q,r)
            if nargin<3
                q=obj.q; r=obj.r;
            end
            v=v-q.*sum(q.*v,3)./r.^2;
        end
        
        % resample v
        function v=resample_aux(obj,q,r)
            if nargin<2
                q=obj.q; r=obj.r;
            end
            v=randn(obj.N,obj.D,obj.w);
            for i=1:obj.w-1
                v(:,i,1:obj.w-i)=0;
            end
            v=obj.mult_proj(v,q,r);
        end
        
        % 1 step of leapfrog
        function [q,v,U,dU]=leapfrog(obj,q,v,dU,h,r,calcU)
            if nargin<5
                h=obj.h; r=obj.r; calcU=true;
            elseif nargin<6
                r=obj.r; calcU=true;
            elseif nargin<7
                calcU=true;
            end
            
            % natural gradient
            ng=obj.mult_proj(dU,q,r);
            
            % forward half step of velocity
            v=v -h/2.*ng;
            
            % full step evolution on sphere along great circles
            v_nom=sqrt(sum(v.^2,3))./r;
            q0 = q;
            cosvt = cos(h.*v_nom); sinvt = sin(h.*v_nom);
            q = q0.*cosvt + v./[ones(obj.N,1),v_nom(:,2:end)].*sinvt;
            v = -q0.*v_nom.*sinvt + v.*cosvt;
%             gcirc = (q+1i*v./[ones(obj.N,1),v_nom(:,2:end)]).*exp(-1i*v_nom.*h);
%             q = real(gcirc); v = imag(gcirc).*v_nom;
            
            % update geometry
            [U,dU] = obj.geom(q,[1,double(~calcU)]);
            ng=obj.mult_proj(dU,q,r);
            
            % backward half step of velocity
            v=v -h/2.*ng;
            
            % adjust the state when necessary
            if any(abs(sum(q(:,:).*v(:,:),2))./obj.D>1e-6)
                [q,v]=obj.adj2S(q,v,r);
            end
        end
        
        % adjust the position to stay on the sphere
        function [q,v]=adj2S(obj,q,v,r)
            if nargin<4
                r=obj.r;
            end
            q=q./sqrt(sum(q.^2,3));
            v=v-q.*sum(q.*v,3);
            q=q.*r;
        end
        
        % standard multiple (wband) spherical HMC
        function [obj,acpt_ind]=std_wbd_SphHMC(obj,rndL,h,L,r)
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
            E_cur=obj.U+(v(:)'*v(:))./2;
            
            % randomize the number of integration steps if asked
            if rndL
                L=ceil(rand*L);
            end
            
            % leapfrogs
            for l=1:L
                [q,v,U,dU]=obj.leapfrog(q,v,dU,h,r,l==L);
            end
            
            % proposed energy
            E_prp=U+(v(:)'*v(:))./2;
            
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
        
        % 1 step of leapfrog for robust algorithm
        function [q,v,U,ng,pow]=leapfrog_rbt(obj,q,v,ng,pow,h,r,endL)
            if nargin<6
                h=obj.h; r=obj.r; endL=true;
            elseif nargin<7
                r=obj.r; endL=true;
            elseif nargin<8
                endL=true;
            end
            
            % forward half step of velocity
            v=v -h/2.*ng;
            
            % full step evolution on sphere along great circles
            v_nom=sqrt(sum(v.^2,3))./r;
            q0 = q;
            cosvt = cos(h.*v_nom); sinvt = sin(h.*v_nom);
            q = q0.*cosvt + v./[ones(obj.N,1),v_nom(:,2:end)].*sinvt;
            v = -q0.*v_nom.*sinvt + v.*cosvt;
%             gcirc = (q+1i*v./[ones(obj.N,1),v_nom(:,2:end)]).*exp(-1i*v_nom.*h);
%             q = real(gcirc); v = imag(gcirc).*v_nom;
            
            % update geometry
            [U,dU] = obj.geom(q,[1,double(~endL)]);
            ng=obj.mult_proj(dU,q,r);
            
            % backward half step of velocity
            v=v -h/2.*ng;
            
            % accumulate the power of force
            if ~endL
                pow=pow + h*(v(:)'*dU(:));
            end
            
            % adjust the state when necessary
            if any(abs(sum(q(:,:).*v(:,:),2))./obj.D>1e-6)
                [q,v]=obj.adj2S(q,v,r);
            end
        end
        
        % robuts multiple (wband) spherical HMC
        function [obj,acpt_ind]=rbt_wbd_SphHMC(obj,rndL,h,L,r)
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
            
            % natural gradient
            ng=obj.mult_proj(dU,q,r);
            
            % accumulate the power of force
            pow=h/2*(v(:)'*dU(:));
            
            % current energy
            E_cur=obj.U - h^2/8.*(dU(:)'*ng(:));
            
            % randomize the number of integration steps if asked
            if rndL
                L=ceil(rand*L);
            end
            
            % leapfrogs
            for l=1:L
                [q,v,U,ng,pow]=obj.leapfrog_rbt(q,v,ng,pow,h,r,l==L);
            end
            
            % accumulate the power of force
            pow=pow +h/2*(v(:)'*ng(:));
            
            % proposed energy
            E_prp=U - h^2/8.*(ng(:)'*ng(:));
            
            % log of Metropolis ratio
            logAP = -E_prp + E_cur + pow;
            % accept or reject the proposal at end of trajectory
            if isfinite(logAP) && (log(rand) < min([0,logAP]))
                obj.q=q; obj.U=U; obj.dU=ng;
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
            E_cur=U0+(v0(:)'*v0(:))./2;
            h=1;
            [q,v,U,dU]=obj.leapfrog(q0,v0,dU0,h);
            E_prp=U+(v(:)'*v(:))./2;
            logAP = -E_prp + E_cur;
            a=2*(exp(logAP)>0.5)-1;
            while a*logAP>-a*log(2)
                h=h*2^a;
                [q,v,U,dU]=obj.leapfrog(q0,v0,dU0,h);
                E_prp=U+(v(:)'*v(:))./2;
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
        
        % adaptive multiple (wband) Spherical HMC
        function [obj,acpt_ind]=adp_wbd_SphHMC(obj,iter,Nadpt)
            % initialization
            q=obj.q; dU=obj.dU;
            v=obj.resample_aux();
            
            % current energy;
            E_cur=obj.U+(v(:)'*v(:))./2;
            
            % leapfrogs
            for l=1:obj.L
%             for l=1:ceil(1/obj.h_adpt.h)
                [q,v,U,dU]=obj.leapfrog(q,v,dU,obj.h_adpt.h);
                % two orthants' rule
                if any(sum(q(:,:).*obj.q(:,:),2)<0)
%                 % symmetric stochastic rule
%                 if any(rand <= (1-mean(q(:,:).*obj.q(:,:),2)./obj.r^2)./2)
                    fprintf('It takes %d steps to exit!\n',l);
                    break;
                end
            end
            
            % proposed energy
            E_prp=U+(v(:)'*v(:))./2;
            
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
            samp=zeros(Nsamp,obj.N,obj.D,obj.w);
            engy=zeros(Niter,1);
            acpt=0; % overall acceptance
            accp=0; % online acceptance
            % run MCMC to collect samples
            fprintf('Running %s wband Spherical HMC...\n',obj.alg);
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
                        [obj,acpt_ind]=obj.std_wbd_SphHMC();
                    case 'robust'
                        [obj,acpt_ind]=obj.rbt_wbd_SphHMC();
                    case 'adaptive'
                        [obj,acpt_ind]=obj.adp_wbd_SphHMC(iter,NBurnIn);
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
                        samp(ceil((iter-NBurnIn)/thin),:,:,:) = obj.q;
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
                samp_diag=samp(:,:,:,end);
                wt=sum(log(abs(samp_diag(:,:))),2);
                wt=exp(wt-median(wt));
                resamp_idx=datasample(1:Nsamp,Nsamp,'replace',true,'weights',wt);
            else
                wt=ones(Nsamp,1);
                resamp_idx=1:Nsamp;
            end
            % save
            if SAVE
                time_lbl=regexprep(num2str(fix(clock)),'    ','_');
                f_name=[obj.alg,'_wbd_SphHMC_dim',num2str(obj.D),'_bdwth',num2str(obj.w),'_',time_lbl];
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