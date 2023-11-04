% Copyright 2021 Fang Bai <fang dot bai at yahoo dot com>
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.




classdef LandmarkSLAMSimulator3D < handle


    properties (Access = public)
        
        GT_Poses = {}
        
        GT_Landmarks = []
 
    end


    properties (Access = public)
              
        dim = 3;
        
        numPoses = 0
        numLandmarks = 0
        
        sigma_sensor_noise = 0.01;
        sigma_rot_noise = 0.01;
        sigma_trans_noise = 0.01;
    end

    
    properties (Access = public)
        
        sensorMinRange = 0
        
        sensorMaxRange = 4

        sensorFOV = (120/180)*pi
                
    end





    methods (Access = public)


        function obj = LandmarkSLAMSimulator3D()

            points_radius = 5;

            num_poses = 20;
            sensor_radius = num_poses / (2*pi);
            N = num_poses/2;
            
            % simulate poses and landmarks
            obj.simulateLandmarks3D(points_radius);
            obj.simulateSensorTrajectory3D(sensor_radius, N);

            obj.numPoses = length(obj.GT_Poses);
            obj.numLandmarks = length(obj.GT_Landmarks);
        end
        

        function [DataCell, VisVecCell, OdometryCell] = generateData (this, sigma_sensor_noise, sigma_rot_noise, sigma_trans_noise)

            if exist("sigma_sensor_noise", "var") && (~isempty(sigma_sensor_noise))
                this.sigma_sensor_noise = sigma_sensor_noise;
            end
            if exist("sigma_rot_noise", "var") && (~isempty(sigma_rot_noise))
                this.sigma_rot_noise = sigma_rot_noise;
            end            
            if exist("sigma_trans_noise", "var") && (~isempty(sigma_trans_noise))
                this.sigma_trans_noise = sigma_trans_noise;
            end

            
            DataCell = cell(1, this.numPoses);
            VisVecCell = cell(1, this.numPoses);
            
            for ii = 1 : this.numPoses
                
                pose = this.GT_Poses{ii};
                R = pose(1:this.dim, 1:this.dim);
                t = pose(1:this.dim, this.dim+1);

                % simulate data
                data = R' * (this.GT_Landmarks - t) + (this.sigma_sensor_noise * randn(this.dim, this.numLandmarks));
                % simulate visibility
                vis = true(1, this.numLandmarks);
                if (1)
                    % points in range [min_range, max_range]
                    vis_min_range = (data(3, :) > this.sensorMinRange);
                    vis_max_range = (data(3, :) < this.sensorMaxRange);
                    % points in fov: angle with respect to zaxis
                    vis_fov = (acos(data(3,:) ./ sqrt(sum(data .* data, 1))) < (this.sensorFOV/2));
                    % visibility
                    vis = logical (vis_min_range .* vis_max_range .* vis_fov);
                end

                DataCell{ii} = data .* vis;
                VisVecCell{ii} = logical(vis);
                
            end

            % simulate odometry
            for ii = 1 : this.numPoses-1
                
                pose1 = this.GT_Poses{ii};
                R1 = pose1(1:this.dim, 1:this.dim);
                t1 = pose1(1:this.dim, this.dim+1);

                pose2 = this.GT_Poses{ii+1};
                R2 = pose2(1:this.dim, 1:this.dim);
                t2 = pose2(1:this.dim, this.dim+1);

                rlv_R = R1' * R2 * SO3.Exp(this.sigma_rot_noise * randn(3, 1));
                rlv_t = R1' * (t2 - t1) + (this.sigma_trans_noise * randn(3, 1));

                OdometryCell{ii} = [rlv_R, rlv_t];
            end
            
        end


    end




    % --- methods to be adapted
    methods (Access = private)

 
        function simulateLandmarks3D (this, points_radius)
        % simulate trajectory
        % output: this.GT_Landmarks
            [X, Y, Z] = sphere(points_radius);

            pole1 = [X(1, 1); Y(1, 1); Z(1, 1)];
            pole2 = [X(end, end); Y(end,end); Z(end, end)];

            X = X(2:end-1, 1:end-1)';
            Y = Y(2:end-1, 1:end-1)';
            Z = Z(2:end-1, 1:end-1)';

            interior = [X(:), Y(:), Z(:)]';

            this.GT_Landmarks = [pole1, interior, pole2];

        end


        function simulateSensorTrajectory3D (this, sensor_radius, N)
        % simulate trajectory
        % output: this.GT_Poses
            [X, Y, Z] = cylinder(sensor_radius, N);

            X = X(:, 1:end-1);
            Y = Y(:, 1:end-1);
            Z = Z(:, 1:end-1);

            X = [X(1, :), X(2, end:-1:1)];
            Y = [Y(1, :), Y(2, end:-1:1)];
            Z = [Z(1, :), Z(2, end:-1:1)];

            centroid = mean(this.GT_Landmarks, 2);

            for ii = 1 : length(X)
                t = [X(ii); Y(ii); Z(ii)];                
                R = this.sensorRotationByPrincipleAxis(centroid - t);
                this.GT_Poses{ii} = [R, t];
            end
        end


    end





    % --- methods (help function)
    methods (Access = private, Static = true)

        function R = sensorRotationByPrincipleAxis (principleAxis)
            
            [U, ~, ~ ] = svd (principleAxis);
            
            R = U(:, [3,2, 1]);
            
            if (R(:, 3)' * principleAxis < 0)
                R(:, 3) = - R(:, 3);
            end
            
            if (det(R) < 0)
                R(:, 1) = - R(:, 1);
            end
            
        end

    end





end