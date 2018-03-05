%%%% ---- Spherical Optimizer ---- %%%%
%%%% Shiwei Lan                    %%%%
%%%% CMS, CalTech                  %%%%
%%%% slan@caltech.edu;             %%%%
%%%% Copyright @ 2017              %%%%
%%%% ---- Spherical Optimizer ---- %%%%

classdef wbandt_SphOpt
    % Optimizer on product of spheres with increasing dimensions.
    % -- modified for applications of processes (time series) --
    % -- speed up by exploring k-band structure of covariance --
    properties
        r=1; % radius of the sphere
        q; N; D; w; % current state and its series length, dimension and diagonal bandwidth
        U; DU; % current potential and its natural (covariant) gradient
        geom; % the function to calculate geometry (potential and gradient)
        h=1; L=1; % the step size and the number of leapfrogs (standard)
        s=0.5; % learning rate in backtrack algorithm
        c=0.1; % control parameter
        print=true; % print the progress
        record=false; % record the progress
    end
    methods
        % constructor
        function obj=wbandt_SphOpt(q,geom,h,L,s,c,PRINT,RECORD)
            obj.q=q; [obj.N,obj.D,obj.w]=size(q);
            obj.r=ones(1,obj.D);
            obj.geom=geom;
            [obj.U,dU]=obj.geom(q,[0,1]);
            obj.DU=obj.mult_proj(dU,obj.q,obj.r); % natural gradient
            if exist('h','var')
                obj.h=h;
            end
            if exist('L','var')
                obj.L=L;
            end
            if exist('s','var')
                obj.s=s;
            end
            if exist('c','var')
                obj.c=c;
            end
            if exist('PRINT','var')
                obj.print=PRINT;
            end
            if exist('RECORD','var')
                obj.record=RECORD;
            end
        end
        
        % projection a vector v onto each of (d-1) sphere S^(i-1)
        function v=mult_proj(obj,v,q,r)
            if nargin<3
                q=obj.q; r=obj.r;
            elseif nargin<4
                r=obj.r;
            end
            v=v-q.*sum(q.*v,3)./r.^2;
        end
        
        % normalize v
        function [v_star,v_nom]=mult_nmlz(obj,v,r)
            if nargin<3
                r=obj.r;
            end
            v_nom=sqrt(sum(v.^2,3))./r;
            v_star=v./[ones(obj.N,1),v_nom(:,2:end)];
        end
        
        % 1 step of leapfrog
        function [q,v,U,DU]=leapfrog(obj,q,v,DU,h,r,calcU)
            if nargin<5
                h=obj.h; r=obj.r; calcU=true;
            elseif nargin<6
                r=obj.r; calcU=true;
            elseif nargin<7
                calcU=true;
            end
            
            % forward half step of velocity
            v=v -h/2.*DU;
            
            % full step evolution on sphere along great circles
            [v_star,v_nom]=obj.mult_nmlz(v,r);
            q0 = q;
            cosvt = cos(h.*v_nom); sinvt = sin(h.*v_nom);
            q = q0.*cosvt + v_star.*sinvt;
            v = -q0.*v_nom.*sinvt + v.*cosvt;
%             gcirc = (q+1i*v./[ones(obj.N,1),v_nom(:,2:end)]).*exp(-1i*v_nom.*h);
%             q = real(gcirc); v = imag(gcirc).*v_nom;
            
            % update geometry
            [U,dU] = obj.geom(q,[1,double(~calcU)]);
            DU=obj.mult_proj(dU,q,r);
            
            % backward half step of velocity
            v=v -h/2.*DU;
            
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
        
        % backtrack
        function obj=backtrack(obj,max_iter,rndL,h,L,r)
            if nargin<2
                max_iter=10; rndL=false; h=obj.h; L=obj.L; r=obj.r;
            elseif nargin<3
                rndL=false; h=obj.h; L=obj.L; r=obj.r;
            elseif nargin<4
                h=obj.h; L=obj.L; r=obj.r;
            elseif nargin<5
                L=obj.L; r=obj.r;
            elseif nargin<6
                r=obj.r;
            end
            
            % backtrack
            for i=1:max_iter
                % initialization
                q=obj.q; DU=obj.DU;
                v=zeros(obj.N,obj.D,obj.w);

                % randomize the number of integration steps if asked
                if rndL
                    L=ceil(rand*L);
                end

                % leapfrogs
                for l=1:L
                    [q,v,U,DU]=obj.leapfrog(q,v,DU,h,r,l==L);
                end
                
                % accept the move if the Armijo-Goldstein is satisfied
                [DU_star,DU_nom]=obj.mult_nmlz(obj.DU);
                if U<obj.U-h.*obj.c* sum(DU_nom(:)) %sum=(DU_star(:)'*obj.DU(:))
                    obj.q=q; obj.U=U; obj.DU=DU;
                    break;
                else
                    h=h.*obj.s;
                end
            end
            
        end
        
        % optimize function
        function [obj,time,f_name]=optimize(obj,Nmax,thld)
            f_name=[];
            if nargin<2
                Nmax=100; thld=1e-3;
            elseif nargin<3
                thld=1e-3;
            end
            % keep track
            minim=zeros(Nmax,obj.N,obj.D,obj.w);
            engy=zeros(Nmax,1);
            % run optimizer
            fprintf('Running wband Spherical Optimizer...\n');
            prog=0.05:0.05:1;
            tic;
            for iter=1:Nmax

                % display the progress
                if obj.print && ismember(iter,floor(Nmax.*prog))
                    fprintf('\t%.0f%% of max iterations completed.\n',100*iter/Nmax);
                end

                % Backtrack to get minimum
                obj=obj.backtrack();

                % keep record
                engy(iter) = obj.U;
                minim(iter,:,:,:) = obj.q;
                if obj.print
                    fprintf('- Objective function value: %.4f\n', engy(iter));
                end
                
                % break if condition satisfied
                if iter>1
                    dif_q=minim(iter,:,:,:)-minim(iter-1,:,:,:);
                    if abs(engy(iter)-engy(iter-1))<thld || max(abs(dif_q(:)))<thld
                        fprintf('- wband Spherical optimization breaks at iteration %d.\n',iter);
                        break;
                    end
                end
                
            end
            % count time
            time=toc;
            fprintf('wband Spherical optimization terminates after %.2f seconds.\n',time);
            % record results
            if obj.record
                time_lbl=regexprep(num2str(fix(clock)),'    ','_');
                f_name=['wbd_SphOpt_dim',num2str(obj.D),'_',time_lbl];
                if exist('result','dir')~=7
                    mkdir('result');
                end
                h=obj.h; L=obj.L;
                save(['result/',f_name,'.mat'],'Nmax','thld','h','L','minim','engy','time');
            end
        end
    end
end