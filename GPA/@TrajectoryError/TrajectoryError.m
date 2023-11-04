classdef TrajectoryError

    methods (Static = true)
        
        function [RMSE_rot, RMSE_pos] = GetAbsoluteTrajectoryErrorByGauge (poseArray1, PoseArrayGauge,  GaugeT)
            
            if ~exist('GaugeT', 'var')
                GaugeT = eye(4);
            end
            
            nsize = length(poseArray1);
            
            pose = PoseArrayGauge{1};
            if (size(pose, 1) == size(pose, 2))
                ndim = size(pose, 1) - 1;
            end
            if (size(pose, 1) + 1 == size(pose, 2))
                ndim = size(pose, 1);
            end
            
            GaugeRot = GaugeT(1:ndim, 1:ndim);
            GaugePos = GaugeT(1:ndim, 1+ndim);
            
            
            rotation_diff = [];
            position_diff = [];
            
            for ii = 1 : length(poseArray1)
                pose1 = poseArray1{ii};
                R = pose1(1:ndim, 1:ndim);
                t = pose1(1:ndim, 1+ndim);
                
                rotation = GaugeRot*R;
                position = GaugeRot*t + GaugePos;
                
                poseGauge = PoseArrayGauge{ii};
                
                gauge_rotation = poseGauge(1:ndim, 1:ndim);
                gauge_position = poseGauge(1:ndim, 1+ndim);
                
                if(ndim == 3)
                    rotation_diff = [rotation_diff,  SO3.Log(rotation' * gauge_rotation)];
                end
                if (ndim == 2)
                    dR = rotation' * gauge_rotation;
                    rotation_diff = [rotation_diff, atan2(dR(2,1), dR(1,1))];
                end
                position_diff = [position_diff,  position-gauge_position];
                
            end
            
            
            RMSE_rot = norm(rotation_diff,  'fro');
            RMSE_pos = norm(position_diff,  'fro');
            
            RMSE_rot = sqrt(RMSE_rot * RMSE_rot / nsize);
            RMSE_pos = sqrt(RMSE_pos * RMSE_pos / nsize);
            
        end
        
        
        
        
        
        function [RMSE_rot, RMSE_pos] = AbsoluteTrajectoryErrorByFirst (poseArray1, PoseArrayGauge)
                        
            nsize = length(poseArray1);
            
            pose = PoseArrayGauge{1};
            if (size(pose, 1) == size(pose, 2))
                ndim = size(pose, 1) - 1;
            end
            if (size(pose, 1) + 1 == size(pose, 2))
                ndim = size(pose, 1);
            end
                        
            tmpPose1 = poseArray1{1};
            tmpPose2 = PoseArrayGauge{2};
            
            T1 = [tmpPose1(1:ndim, 1:ndim), tmpPose1(1:ndim, 1+ndim);  [zeros(1, ndim), 1]];
            T2 = [tmpPose2(1:ndim, 1:ndim), tmpPose2(1:ndim, 1+ndim);  [zeros(1, ndim), 1]];
            
            GaugeT =  T2 * inv(T1);
            
            
            GaugeRot = GaugeT(1:ndim, 1:ndim);
            GaugePos = GaugeT(1:ndim, 1+ndim);
            
            
            rotation_diff = [];
            position_diff = [];
            
            for ii = 1 : length(poseArray1)
                pose1 = poseArray1{ii};
                R = pose1(1:ndim, 1:ndim);
                t = pose1(1:ndim, 1+ndim);
                
                rotation = GaugeRot*R;
                position = GaugeRot*t + GaugePos;
                
                poseGauge = PoseArrayGauge{ii};
                
                gauge_rotation = poseGauge(1:ndim, 1:ndim);
                gauge_position = poseGauge(1:ndim, 1+ndim);
                
                if(ndim == 3)
                    rotation_diff = [rotation_diff,  SO3.Log(rotation' * gauge_rotation)];
                end
                if (ndim == 2)
                    dR = rotation' * gauge_rotation;
                    rotation_diff = [rotation_diff, atan2(dR(2,1), dR(1,1))];
                end
                position_diff = [position_diff,  position-gauge_position];
                
            end
            
            
            RMSE_rot = norm(rotation_diff,  'fro');
            RMSE_pos = norm(position_diff,  'fro');
            
            RMSE_rot = sqrt(RMSE_rot * RMSE_rot / nsize);
            RMSE_pos = sqrt(RMSE_pos * RMSE_pos / nsize);
            
        end        
        
        
        
        
        function [RMSE_rot, RMSE_pos] = RelativePoseError (poseArray1, PoseArrayGauge)
            
            nsize = length(poseArray1);
            
            pose = PoseArrayGauge{1};
            if (size(pose, 1) == size(pose, 2))
                ndim = size(pose, 1) - 1;
            end
            if (size(pose, 1) + 1 == size(pose, 2))
                ndim = size(pose, 1);
            end
                        
                        
            rotation_diff = [];
            position_diff = [];
            
            
            cnt = 1;
            
            R0 = poseArray1{cnt}(1:ndim, 1:ndim);
            t0 = poseArray1{cnt}(1:ndim, 1+ndim);
            
            gauge_R0 = PoseArrayGauge{cnt}(1:ndim, 1:ndim);
            gauge_t0 = PoseArrayGauge{cnt}(1:ndim, 1+ndim);
            
            
            for cnt = 2 : nsize
                
                R = poseArray1{cnt}(1:ndim, 1:ndim);
                t = poseArray1{cnt}(1:ndim, 1+ndim);
                               
                gauge_R = PoseArrayGauge{cnt}(1:ndim, 1:ndim);
                gauge_t = PoseArrayGauge{cnt}(1:ndim, 1+ndim);
               
                rotation = R0' * R;
                position = R0' * (t - t0);                
                
                gauge_rotation = gauge_R0' * gauge_R;
                gauge_position = gauge_R0' * (gauge_t - gauge_t0);
                                
                R0 = R;  t0 = t;
                gauge_R0 = gauge_R;  gauge_t0 = gauge_t;
                
                if(ndim == 3)
                    rotation_diff = [rotation_diff,  SO3.Log(rotation' * gauge_rotation)];
                end
                if (ndim == 2)
                    dR = rotation' * gauge_rotation;
                    rotation_diff = [rotation_diff, atan2(dR(2,1), dR(1,1))];
                end
                position_diff = [position_diff,  position-gauge_position];
                
            end
            
            
            RMSE_rot = norm(rotation_diff,  'fro');
            RMSE_pos = norm(position_diff,  'fro');
            
            RMSE_rot = sqrt(RMSE_rot * RMSE_rot / nsize);
            RMSE_pos = sqrt(RMSE_pos * RMSE_pos / nsize);
            
        end                
        
        
    end
end

