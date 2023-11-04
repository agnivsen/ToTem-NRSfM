classdef NonRigidTransformations<handle
    properties(Access=private)
        Linv = [];
        W = [];
        A = [];
        d = [];
        base = [];
        rhoThreshold = 0.0000001;
    end
    
    methods(Access=public)
        %%
        %%%%%%%%%%%%%%%%%%%%
        %  Compute the 3D TPS parameters corresponding to target positions
        %  for the control handles
        %
        % @param: base - the input control handles for the TPS
        % @param: N - number of 3D points in data
        %
        % @return: Object of the TPS3D class
        %%%%%%%%%%%%%%%%%%%%
        function [obj] = NonRigidTransformations(base)
            N = size(base, 1);
            target = base;
            
            if(size(base,2) == 2)
                target = [target ones(N,1)];
                base = [base ones(N,1)];
            end
            
            K = zeros(N, N, 'double');
            P = ones (N, 4, 'double');
            O = zeros(4, 4, 'double');
            Y = zeros((N+4), 4, 'double');
            
            for i = 1:N
                for j = 1:N
                    if(i~=j)
                        x_n = base(i,1) - base(j,1);
                        y_n = base(i,2) - base(j,2);
                        z_n = base(i,3) - base(j,3);
                        K(i,j) = obj.rho(x_n, y_n, z_n);
                    end
                end
                P(i,1) = base(i,1);
                P(i,2) = base(i,2);
                P(i,3) = base(i,3);
                P(i,4) = 1.0;
                Y(i,1) = target(i,1);
                Y(i,2) = target(i,2);
                Y(i,3) = target(i,3);
                Y(i,4) = 1.0;
            end
            
            Kinv = pinv(K);
            xi = Kinv - ( Kinv *  P * pinv(P' * Kinv * P) * P' * Kinv );
            
            KP = horzcat(K,P);
            PTO = horzcat(P',O);
            L = vertcat(KP,PTO);
            
            Linv = pinv(L);
            W = Linv*Y;
            A = W(1:N,:);
            d = W(N+1:end,:);
            
            obj.Linv = Linv; obj.W = W; obj.A = A; obj.d = d;  obj.base = base;          
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%
        %  Update the TPS based on updated positions of the control handle
        %
        % @param: target - the updated target position for the control
        %                                   handles of the TPS
        % @param: Linv - inverse of design matrix, pre-computed
        %%%%%%%%%%%%%%%%%%%%
        function [obj] = update_tps(obj, target)
            N = size(target, 1);
            if(size(target,2) == 2)
                target = [target ones(N,1)];
            end            
            Y = zeros((N+4), 4, 'double');
            
            for i = 1:N
                Y(i,1) = target(i,1);
                Y(i,2) = target(i,2);
                Y(i,3) = target(i,3);
                Y(i,4) = 1.0;
            end
            
            obj.W = obj.Linv*Y;
            obj.A = obj.W(1:N,:);
            obj.d = obj.W(N+1:end,:);
        end
        
        function [interpolated] = interpTPS(obj, interpolationCandidates)
            N = size(interpolationCandidates, 1);
            if(size(interpolationCandidates,2) == 2)
                interpolationCandidates = [interpolationCandidates ones(N,1)];
            end      
            
            cntrlN = size(obj.base,1);
            interpolated = zeros(N,3);
            
            for i = 1:N
                P_t = ones(1, 4);
                P_t(1,1) = interpolationCandidates(i,1);
                P_t(1,2) = interpolationCandidates(i,2);
                P_t(1,3) = interpolationCandidates(i,3);
                Phi_t = zeros(1, cntrlN);

                for j = 1:cntrlN
                    Phi_t(1,j) = obj.rho(  obj.base(j,1) -  P_t(1,1),  ...
                                                obj.base(j,2) -  P_t(1,2), ...
                                                obj.base(j,3) -  P_t(1,3)  );
                end
                fp = (P_t*obj.d) + (Phi_t*obj.A);
                interpolated(i,:) = fp(1:3);
            end            
            
            if(size(interpolationCandidates,2) == 2)
                interpolated(:,3) = [];
            end
        end

        function [interpolated] = interpTPSNormal(obj, interpolationCandidates)
            N = size(interpolationCandidates, 1);
            if(size(interpolationCandidates,2) == 2)
                interpolationCandidates = [interpolationCandidates ones(N,1)];
                warning('Normal interpolation only for R2 -> R3 warps!');
            end      
            
            cntrlN = size(obj.base,1);
            interpolated = zeros(N,3);
            
            for i = 1:N
                P_t = ones(1, 4);
                P_t(1,1) = interpolationCandidates(i,1);
                P_t(1,2) = interpolationCandidates(i,2);
                P_t(1,3) = interpolationCandidates(i,3);
                Phi_t_x = zeros(1, cntrlN);
                Phi_t_y = zeros(1, cntrlN);

                for j = 1:cntrlN
                    Phi_t_j = obj.rho(  obj.base(j,1) -  P_t(1,1),  ...
                                                obj.base(j,2) -  P_t(1,2), ...
                                                obj.base(j,3) -  P_t(1,3)  );
                    Phi_t_x(1,j) = -(obj.base(j,1) -  P_t(1,1))/Phi_t_j;
                    Phi_t_y(1,j) = -(obj.base(j,2) -  P_t(1,2))/Phi_t_j;

                end
                Nx = ([1 0 0 0]*obj.d) + (Phi_t_x*obj.A);
                Ny = ([0 1 0 0]*obj.d) + (Phi_t_y*obj.A);
                Nx = Nx(1:3); Ny = Ny(1:3);
                crossNxNy = cross(Nx,Ny);
                N = crossNxNy./norm(crossNxNy);
                interpolated(i,:) = N;
            end            
        end        
        
    end
    
    methods(Access=private)        
        function Ur = rho(obj, x,y,z)
            Ur = sqrt((x*x) + (y*y) + (z*z));
            if(abs(Ur) < obj.rhoThreshold)
                Ur = 0.0001;
            end
        end
    end
end