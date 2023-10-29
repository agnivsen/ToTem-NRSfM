classdef SphereAlignment < handle
    methods
%%             Fits a sphere to given pointcloud
% % %         @param: Input: [Nx3] pointcloud that needs to be aligned
% % %         Output:  
%%%           @returns: aligned_cloud - returns the cloud aligned with the
%%%           origin of R3
% % %         @returns: r - radius of the fitted circle
% % %         @returns: C - center of fitted circle
%%%           @returns: flag - 0: if SDP failed, 1: if success
        function[ aligned_cloud, r, C, T, flag ] = align_pointcloud(obj, cloud)

            try 
                pars.fid = 0;
                K.l = 2;
                sedumi(ones(1,2), 1, 0, K, pars); 
            catch 
                fprintf('<strong>SeDuMi not installed</strong>\n[Installing SeDuMi now...]\n\n');
                pause(1);
                cd './Dependencies/SeDuMi_1_3/'
                install_sedumi;
                cd './../../'
                fprintf('\n\n<strong>SeDuMi installed</strong>\n');
            end
            mpol cx cy cz rSquared;

            val = 0;
            N_data = size(cloud);
            N_data = N_data(1);

            for index = 1:N_data
                X = cloud(index,1);
                Y = cloud(index,2);
                Z = cloud(index,3);

                d1 = (X - cx).^2 + (Y - cy).^2 + (Z - cz).^2;
                res = (d1 - rSquared).^2;
                val = val + res;
            end
            
            val = val/(N_data.^2);

            P = msdp(min(val));
            [status,obj] = msol(P);
            status;
            
            flag = 1;

%             status = -1;

            if(status <=0) %%% Contingency measure of sphere fitting fails
                
                % % % ---: Our method
                warning('Sphere fitting failed, resorting to our custom heuristics - centering and computing radius!');
                C = mean(cloud);
                centered = cloud - C;
                norm(centered);
                r = -1;
                for index = 1:N_data
                    dist = norm(centered(index,:));
                    if(dist > r)
                        r = dist;
                    end
                end
            else
                x = double([cx cy cz rSquared]);
                C = x(1:3);
                r = sqrt(x(4));
            end
            
            aligned_cloud = zeros(N_data,3);
            
            scale = 1/r;
            
            for index = 1:N_data
                aligned_cloud(index,:) = scale*(cloud(index,:) - C);
            end
            
             
            T = eye(4);
            S = eye(4,4);
            S(1,1) = 1/scale; S(2,2) = 1/scale; S(3,3) = 1/scale;
            T = S*T;
            T(1:3,4) = C.';
            
        end
    end
end
        
        
