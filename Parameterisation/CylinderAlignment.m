classdef CylinderAlignment < handle
    methods
% %         @param: Input: [Nx3] pointcloud that needs to be aligned
% %         Output:  
% %         @returns: aligned_cloud - [Nx3] pointcloud aligned with Y-axis
% %         @returns: R - rotation matrix to rotate data (at frame A) to be aligned to Y axis, i.e. {}^YR_A
% %         @returns: C - center of orthographically projected circle. To map R2 |-> R3, do [cx,cy] |-> [cx, 0, cy]
% %         @returns: T - tranformation from Y-aligned canonical pose to the previous pose of data (which was given as input)
function[aligned_cloud, R, C, T] = align_pointcloud(obj, cloud, referenceVector)

            if exist('referenceVector', 'var')
                direction_vector = referenceVector;
                centroid_initial = mean(cloud);
                cloud_centered = cloud - centroid_initial;
            else
                centroid_initial = mean(cloud);
                cloud_centered = cloud - centroid_initial;
                [~,~,V] = svd(cloud_centered);
                direction_vector = V(:,1);
            end
            
            
            Y_direction = [0;1;0];
            
            angle_with_Y_axis = rad2deg(acos(dot(direction_vector, Y_direction)));
            
            if(abs(angle_with_Y_axis) > 100) % if angle between vector & Y-axis greater tha 100 degrees (opp. aligned), flip direction vector
                direction_vector = -direction_vector;
            end
            
            v = cross(direction_vector, Y_direction);
            
            s = norm(v);
            c = dot(direction_vector, Y_direction);
            v_hat =[0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0 ];
            
            R = eye(3) + v_hat + (v_hat*v_hat*(1-c)/s.^2);
            N_data = size(cloud_centered);
            N_data = N_data(1);
            ortho_cloud = zeros(N_data,3);
            rot_ortho_cloud = zeros(N_data,3);
            aligned_cloud = zeros(N_data,3);
            
            for index = 1:N_data
                q = cloud_centered(index,:).';
                res = dot(q,direction_vector);
                % Orthogonal projection here
                pQ = q - ((dot(q,direction_vector))/(dot(direction_vector,direction_vector))*direction_vector);
                ortho_cloud(index,:) = pQ.';
                pQ = R*pQ;
                rot_ortho_cloud(index,:) = pQ.';
                q = R*q;
                aligned_cloud(index,:) = q.';
            end

            rot_ortho_cloud(:,2) = [];
            
            [C,~] = obj.circle_fit(rot_ortho_cloud);
            
            
            for index = 1:N_data
                 aligned_cloud(index,:) =  aligned_cloud(index,:) - [C(1) 0 C(2)]; %accounting for displacement
            end
            
            range = [max(aligned_cloud(:,2)) min(aligned_cloud(:,2))];
            
            [regParams,~,~]=absor(aligned_cloud.',cloud.','doScale',false);
            T = regParams.M;
            
            scale = getLargestDiagonal(aligned_cloud);

            aligned_cloud = (1/scale).*aligned_cloud;
            t_temp = T(1:3,4);
            T = scale*T;
            T(1:3,4) = t_temp;

        end
        
        
        function [C,r] = circle_fit(obj, cloud)            
            N_data = size(cloud);
            N_data = N_data(1);
            C = zeros(N_data,3);
            d = zeros(N_data,1);
            
            for index = 1:N_data
                X = cloud(index,1);
                Y = cloud(index,2);                
                C(index,1) = 1;
                C(index,2) = 2*X;
                C(index,3) = 2*Y;
                d(index,1) = X.^2 + Y.^2;
            end
            
            x = lsqlin(C,d);
            rho_ = x(1);
            Cx = x(2);
            Cy = x(3);
            r = sqrt(rho_ + (Cx.^2 + Cy.^2));
            
            C = [Cx Cy];
            
        end
        function plotcircle(obj,r,x,y)
            th = 0:pi/100:2*pi;
            f = r * exp(j*th) + x+j*y;
            plot(real(f), imag(f),'LineWidth',3,'Color','blue');
        end
    end
end


