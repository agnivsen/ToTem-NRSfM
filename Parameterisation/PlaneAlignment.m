classdef PlaneAlignment < handle
    methods
%%             Fits a plane to given pointcloud using SVD
%%%         @param: Input: [Nx3] pointcloud that needs to be aligned
%%%                               Output:  
%%%         @returns: aligned_cloud - returns the cloud aligned with the
%%%                               origin of R3
%%%         @returns: R - Rotation matrix of the fitted plane w.r.t a
%%%                               plane at origin
%%%         @returns: C - center of fitted plane
%%%         @returns: T - transformation matrix from canonical pose to
%%%                               actual pose of fitted plane

        function[ aligned_cloud, R, C, T ] = align_pointcloud_svd(obj, cloud)
            [~,~,V] = svd(cloud);

            first_tangent_vector = V(:,1);
            second_tangent_vector = V(:,2);
            normal = cross(first_tangent_vector,second_tangent_vector);
            
            camera_direction = [0 0 1];
            
            r = vrrotvec(normal,camera_direction);
            R = vrrotvec2mat(r);
            R = eye(3);
            C = mean(cloud);
            T = horzcat(R, C.');
            T = vertcat(T, [0 0 0 1]);
            
            aligned_cloud = R.'*cloud.';
            aligned_cloud = aligned_cloud(1:3,:)  - C.';
            aligned_cloud = aligned_cloud.';
            warning('Need more sophisticated planar fitting algorithm, turned off for now');

        end

%%              Fits a plane to given pointcloud using Linear Least Squares
%%%         @param: Input: [Nx3] pointcloud that needs to be aligned
%%%                               Output:  
%%%         @returns: aligned_cloud - returns the cloud aligned with the
%%%                               origin of R3
%%%         @returns: R - Rotation matrix of the fitted plane w.r.t a
%%%                               plane at origin; origin to fitted plane
%%%                               rotation
%%%         @returns: C - center of fitted plane
%%%         @returns: T - transformation matrix from canonical pose to
%%%                               actual pose of fitted plane

        function[ aligned_cloud, R, C, T ] = align_pointcloud(obj, cloud)
            N_cloud = size(cloud);
            N_cloud = N_cloud(1);
            
            C = cloud;
            d = ones(N_cloud,1);
            normal = lsqlin(C,d);
            normal = normal/norm(normal);
            
            camera_direction = [0 0 1];
            
            dist_sum = 0;

            for index = 1:N_cloud
                dist_sum = dist_sum + ( cloud(index,1)*normal(1) + cloud(index,2)*normal(2) + cloud(index,3)*normal(3));
            end
            
            avgDist = (dist_sum/N_cloud);
            C = normal*avgDist;
            C = C.';
            
            C = mean(cloud);
            
            r = vrrotvec(camera_direction, normal);
            R = vrrotvec2mat(r);
            T = horzcat(R, C.');
            T = vertcat(T, [0 0 0 1]);
            
            aligned_cloud = R.'*cloud.';
            aligned_cloud = aligned_cloud(1:3,:)  - (R.'*C.');
            aligned_cloud = aligned_cloud.';
        end
        
    end
end

