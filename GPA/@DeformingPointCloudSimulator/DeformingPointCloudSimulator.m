classdef DeformingPointCloudSimulator < handle
     
    properties
       
        GT_Poses = {}
        
        GT_InitialLandmarks = []
        
        GT_DeformedLandmarks = {}

        dim = 3;
        
        sigma_deform = 0.1
        
        numPoses = 10
        
        numLandmarks = 10
        
        sigma_noise = 0.01;

    end

    
    properties
        
        sensorMinRange = 0
        
        sensorMaxRange = inf
        
        DefTarget_Landmarks = []
        
    end
    

    properties

        warpModelFlag = 'TPS';

        TPS = struct ( 'IntrisicMatr', [], 'ControlPoints', [], 'InternalSmoothParam', 0, 'Weights', [], 'dimFeatureSpace', [])

    end

    
    
    
    methods (Access = public)
        function obj = DeformingPointCloudSimulator(numPoses, numLandmarks)
            %DEFORMINGPOINTCLOUDSIMULATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.numPoses = numPoses;
            obj.numLandmarks = numLandmarks;
            if exist('./DefGPA', 'dir')
                addpath('./DefGPA');
            end
            if exist('../DefGPA', 'dir')
                addpath('../DefGPA');
            end
        end
        
        
        
        function [DataCell, VisVecCell] =  SimulateMeasurements (this)            
            this.SimulateTrajectoryAndLandmarks();            
            
            this.TPS.dimFeatureSpace = 3^this.dim;            
            this.initTPSWarpModel();

            this.GenerateTargetLandmarks ();
            
            DataCell = cell(1, this.numPoses);
            VisVecCell = cell(1, this.numPoses);
            
            for ii = 1 : this.numPoses
                
                [measurement_data,  vis_logical] = SimulationPoseMeasurement (this, ii);
                
                DataCell{ii} = measurement_data;
                
                VisVecCell{ii} = logical(vis_logical);
                
            end
            
        end
        
        
        
        
    end
    
    
    methods (Access = private)
        

        
        
        
        function [measurements,  vis] = SimulationPoseMeasurement (this, pose_id )
            pose = this.GT_Poses{pose_id};
            rotation = pose(1:this.dim, 1:this.dim);
            position = pose(1:this.dim, 1+this.dim);
            
            vis = this.getVisiblePoints (pose_id);

            % this use linear interpolation to ensure sequential smoothness
            %deformed_landmarks = this.GetDeformedPointCloudByVectorDeformationField(pose_id);
            %deformed_landmarks = this.GetDeformedPointCloudByTPS(pose_id);
    
            % this use independent smooth-deformation for each case
            %deformed_landmarks = this.GetDeformedPointCloudByVectorDeformationField();
            deformed_landmarks = this.GetDeformedPointCloudByTPS();

            this.GT_DeformedLandmarks{pose_id} = deformed_landmarks;

            additive_noise =  randn(size(deformed_landmarks, 1), size(deformed_landmarks, 2));
            additive_noise = min(additive_noise, 3);
            additive_noise = max(additive_noise, -3);

            measurements = rotation' * (deformed_landmarks - position) + this.sigma_noise * additive_noise; 
            measurements = measurements .* vis;
        end
        
        
        function vis = getVisiblePoints (this, pose_id)
            
            % transform points into the camera frame
            pose = this.GT_Poses{pose_id};
            
            rotation = pose(1:this.dim, 1:this.dim);
            position = pose(1:this.dim, 1+this.dim);
     
            ptsLocalFrame = rotation' * (this.GT_InitialLandmarks - position);
            
            distances = sqrt(sum(ptsLocalFrame.*ptsLocalFrame, 1));
            
            % points in range [min_range, max_range]
            vis_min_range = (distances > this.sensorMinRange);
            vis_max_range = (distances < this.sensorMaxRange);
            
            vis = logical (vis_min_range .* vis_max_range);
        end        
        

        function output_point_cloud = GetDeformedPointCloudByVectorDeformationField (this, pose_id)
            if exist("pose_id", 'var')
                alpha = pose_id / this.numPoses;
                output_point_cloud = (1-alpha) * this.GT_InitialLandmarks(:, logical(vis)) + alpha * this.DefTarget_Landmarks(:, logical(vis));
            else
                smooth_deformation_field =  this.SimulateSmoothDeformationField (this.GT_InitialLandmarks);
                output_point_cloud= this.GT_InitialLandmarks + this.sigma_deform * smooth_deformation_field;
            end
        end


        function output_point_cloud = GetDeformedPointCloudByTPS (this, pose_id)
            if exist("pose_id", 'var') && exist("vis", 'var')
                alpha = pose_id / this.numPoses;
                output_point_cloud = (1-alpha) * this.GT_InitialLandmarks(:, logical(vis)) + alpha * this.DefTarget_Landmarks(:, logical(vis));
            else
                output_point_cloud = zeros(size(this.GT_InitialLandmarks, 1), size(this.GT_InitialLandmarks, 2));
                new_control_points_positions = this.TPS.ControlPoints + this.sigma_deform * max(min(randn(size(this.TPS.ControlPoints, 1),  size(this.TPS.ControlPoints, 2)), 3), -3);
                for ii = 1 : size(this.GT_InitialLandmarks, 2)
                    nlv = ThinPlateSpline.getEtaTPS (this.TPS.ControlPoints, this.GT_InitialLandmarks(:, ii));
                    output_point_cloud(:, ii) = new_control_points_positions * this.TPS.IntrisicMatr' * nlv;
                end
            end
        end

        
        
    end
    
    
    
    methods (Access = private)
        
        function SimulateTrajectoryAndLandmarks(this)
            
            sensor_radius = this.numPoses / (2*pi);            
            pts_radius = 0.75 * sensor_radius;            
            zscale = 0.5 * pts_radius;
            
            this.sensorMinRange = 1.0 * norm([sensor_radius-pts_radius,  zscale]);
            this.sensorMaxRange = 1.0 * norm([sensor_radius,  pts_radius-zscale]);

            this.GT_InitialLandmarks = this.SimuPointsBySphere (this.numLandmarks, pts_radius, [0, 0, 0]);            

            % simulate sensor positions
            [CX, CY, CZ] = cylinder(sensor_radius, round(this.numPoses)+0);
            CZ = zscale * CZ;
            pose_positions = [CX(1, :);
                                             CY(1, :);
                                             CZ(2, :)];
            if (this.dim == 2)
                pose_positions = pose_positions ([1,2], :);
            end
                                                
            % simulate sensor orientations
            for ii = 1 : this.numPoses
                
                xaxis = pose_positions(:, ii+1) - pose_positions(:, ii);
                xaxis = xaxis / norm(xaxis);
                
                if (this.dim == 3)
                    zaxis = [0; 0; 1];
                    yaxis = cross(zaxis, xaxis);
                    rotation = [xaxis, yaxis, zaxis];
                end
                
                if (this.dim == 2)
                    yaxis = [-xaxis(2); xaxis(1)];
                    rotation = [xaxis,  yaxis];
                end

                position = pose_positions(:, ii);

                this.GT_Poses{ii} = [rotation,   position];
            end
        end


        function initTPSWarpModel (this)
            this.TPS.ControlPoints =  ThinPlateSpline.setControlPoints (this.GT_InitialLandmarks, this.TPS.dimFeatureSpace);
            this.TPS.InternalSmoothParam = 0; % diagonal elements of TPS kernel matrix. default = 0
            this.TPS.IntrisicMatr = ThinPlateSpline.setTPSIntrisicMatr (this.TPS.ControlPoints, this.TPS.InternalSmoothParam);
        end
        

        function GenerateTargetLandmarks (this)
            smooth_deformation_field =  this.SimulateSmoothDeformationField (this.GT_InitialLandmarks);
            this.DefTarget_Landmarks = this.GT_InitialLandmarks + this.sigma_deform * smooth_deformation_field;
        end
        
    end
    
    
    methods (Static = true, Access = private)
        
        function pts = SimuPointsBySphere (NumPoints, radius, center)
            P = zeros(NumPoints,3) ;
            for i=1:NumPoints
                x1 = rand;
                x2 = rand;
                Px = sqrt(x1*(2-x1))* cos(2*pi*x2);
                Py = sqrt(x1*(2-x1))* sin(2*pi*x2);
                Pz = 1-x1;
                Px = Px * radius + center(1);
                Py = Py * radius + center(2);
                Pz = Pz * radius + center(3);
                P(i,:) = [Px, Py, Pz];
            end
            pts = P';
        end



        function smooth_deformation_field =  SimulateSmoothDeformationField (QueryPoints, n_cut_off, nSampleSize)
            if ~exist("QueryPoints", 'var') || (exist("QueryPoints", 'var') && isempty(QueryPoints))
                dim = 3;
            else
                dim = size(QueryPoints, 1);
            end
            if ~exist('n_cnt', 'var')
                n_cut_off = 4;
            end
            if ~exist('nSampleSize', 'var')
                nSampleSize = 100;
            end

            % Random numbers by Gaussian Distribution
            %X = randn(dim, nSampleSize); X = max(X, -3); X = min(X, 3);
            % Random numbers by Uniform Distribution            
            X = 2*(rand(dim, nSampleSize) - 0.5); 
            % Fourier Transform
            F = fft(X, [], 2);
            % Cut-off_frequency to filter out high frequency variations
            F ( :,  (1+n_cut_off) : (end-n_cut_off+1) ) = 0;
            % Inverse Fourier Transform to obtain Smooth Random Numbers
            SmoothX = ifft (F, [], 2, 'symmetric');

            if exist('QueryPoints', 'var')
                % USE INTERPOLATION FOR QUERY POINTS
                RangeMinMax = [min(QueryPoints, [], 2),  max(QueryPoints, [], 2)];
                lengthXYZ = RangeMinMax (:, 2) - RangeMinMax(:, 1);
                scale = max (lengthXYZ) / (nSampleSize - 1);
                % It is important to USE THE SAME SCALE to ensure uniform deformation in X-Y-Z coordinates
                int_vec = 0 : (nSampleSize-1);
                int_vec(1) = int_vec(1) - 0.1; 
                int_vec(end) = int_vec(end) + 0.1;
                vals = scale * ones(dim ,1) * int_vec + RangeMinMax(:, 1);
                smooth_deformation_field = zeros( size(QueryPoints, 1), size(QueryPoints, 2) );

                % obtain the deformation at query-points by interpolation
                for dd = 1 : size(QueryPoints, 1)
                    smooth_deformation_field(dd, :) = interp1(vals(dd, :), SmoothX(dd, :), QueryPoints(dd, :));
                end
            else
                smooth_deformation_field = SmoothX;
            end

            if nnz (isnan (smooth_deformation_field))
                fprintf(2, "find error\n");
                smooth_deformation_field = [];
            end

            nSampleSize = size(QueryPoints, 2);
            X = 2*(rand(dim, nSampleSize) - 0.5); 
            smooth_deformation_field = X;

        end


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

