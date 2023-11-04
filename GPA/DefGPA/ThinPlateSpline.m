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


classdef ThinPlateSpline
% 
% Thin-Plate-Spline (TPS) implementation as a static class
% maintainer: Fang Bai <fang.bai@yahoo.com>
% created on April 23, 2021, Clermont-Ferrand
% 
% Methods:
%       dist = TPS_BasisFunc ( dp, dim )
%       controlPoints = setControlPoints (Data, dimFeatureSpace)
%       IntrisicMatr = setTPSIntrisicMatr (ControlPoints, InternalSmoothParam)
%       nlv = getEtaTPS (ControlPoints, pVal)
    
    methods (Access = public, Static = true)

        % Thin-Plate-Spline basis function for 2D and 3D case
        % it's a radial distance function
        % dp = the differnece of two points
        function dist = TPS_BasisFunc ( dp, dim )
            
            if (dim == 3)
                
                r = norm (dp);
                
                dist = - r;
                
            end
            
            if (dim == 2)
                
                rr = dp(1) * dp(1) + dp(2) * dp(2) + 1e-20;
                
                dist = rr * log (rr);
                
            end
            
        end
        
        
        % a strategy to choose TPS control points/centers
        function controlPoints = setControlPoints (Data, dimFeatureSpace)
            
            mu = mean(Data,2);

            [U, ~, ~] = svd(Data - mu, 'econ');

            ProjData = U' * (Data - mu);
            
            min_xyz = min (ProjData, [], 2);
            max_xyz = max (ProjData, [], 2);
            spn_xyz = max_xyz - min_xyz;
            dim = size(Data, 1);
            
            numPointsPerAxis= nthroot(dimFeatureSpace, dim);
            numPointsPerAxis2= nthroot(dimFeatureSpace/2, dim-1);
           
            UseRandomAssignment = false;
            
            if (floor(numPointsPerAxis) == numPointsPerAxis)
                
               numPointsPerAxisX = numPointsPerAxis;
                numPointsPerAxisY = numPointsPerAxis;
                numPointsPerAxisZ = numPointsPerAxis;
                            
            elseif dim==3 && (floor(numPointsPerAxis2) == numPointsPerAxis2)
                
                numPointsPerAxisX = numPointsPerAxis2;
                numPointsPerAxisY = numPointsPerAxis2;
                numPointsPerAxisZ = 2;
                
            else
                
                UseRandomAssignment = true;
                
            end
                  
            
            if ~UseRandomAssignment

%                 numPointsPerAxisX = 5;
%                 numPointsPerAxisY = 5;
%                 numPointsPerAxisZ = 1;
                
                xaxis = kron(1:numPointsPerAxisX, ones(1, numPointsPerAxisY)) ./ (numPointsPerAxisX+1);
                yaxis = kron(ones(1, numPointsPerAxisX), 1:numPointsPerAxisY) ./ (numPointsPerAxisY+1);
                
                selectMatr = [xaxis; yaxis];
                
                if (dim == 3)
                    
                    selectMatr = kron(ones(1,numPointsPerAxisZ), selectMatr);
                    
                    zaxis = kron(1:numPointsPerAxisZ, ones(1, numPointsPerAxisX * numPointsPerAxisY)) ./ (numPointsPerAxisZ+1);
                    
                    selectMatr = [selectMatr; zaxis];
                    
                end

            else
                
                selectMatr = rand(dim, dimFeatureSpace);
                
            end
            
            controlPointsProjCoords = min_xyz + (selectMatr .* spn_xyz);
            
            controlPoints = U * controlPointsProjCoords + mu;
            
            if (0)
                figure('Name', 'TPS control-points assignment', 'Position', [1000, 0, 600, 500]);
                clf;
                hold on;
                if (dim == 3)
                    scatter3(Data(1,:), Data(2,:), Data(3,:));
                    scatter3(controlPoints(1,:), controlPoints(2,:), controlPoints(3,:), 'g*');
                    view(3);
                end
                if (dim == 2)
                    scatter(Data(1,:), Data(2,:));
                    scatter(controlPoints(1,:), controlPoints(2,:), 'g*');
                end                
                hold off;
            end

        end
        
        
        
        % copmute the TPS intrinsic matrix
        function IntrisicMatr = setTPSIntrisicMatr (ControlPoints, InternalSmoothParam)
            
            [dim, dimFeatureSpace] = size(ControlPoints);
            
            K = zeros (dimFeatureSpace, dimFeatureSpace);
            
            for it_r = 1 : dimFeatureSpace
                
                for it_c = it_r+1 : dimFeatureSpace
                    
                    dp = ControlPoints(:, it_r) - ControlPoints(:, it_c);
                    
                    K (it_r, it_c) = ThinPlateSpline.TPS_BasisFunc( dp , dim);
                    
                end
                
            end
            
            K = K + K';
            
            for it_d = 1 : dimFeatureSpace
                
                K (it_d, it_d) = InternalSmoothParam;
                
            end
            
            InvK = pinv (K);
            
            C = [ ControlPoints; ones(1, dimFeatureSpace) ];
            
            H = pinv(C * InvK * C') * C * InvK;  % in case control point-cloud is flat
            
            IntrisicMatr = [ InvK - InvK * C' * H;  H ];
            
        end
        
        
        
        function nlv = getEtaTPS (ControlPoints, pVal)
            
            [dim, dimFeatureSpace] = size(ControlPoints);
            
            nlv = ones (dimFeatureSpace + dim + 1, 1);
            
            dmatr = ControlPoints - pVal * ones (1, dimFeatureSpace);
            
            for k = 1 : dimFeatureSpace
                
                nlv (k) = ThinPlateSpline.TPS_BasisFunc (dmatr(:, k), dim);
                
            end
            
            nlv (dimFeatureSpace+1 : dimFeatureSpace+dim) = pVal;
            
        end
        
        

    end

end