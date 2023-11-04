

clear all
close all
clc


DIM = 2; % choose dimenion, 2 for 2D;  3 for 3D

dimFeatureSpace = 10;
InternalSmoothParam = 0;
SmoothParam = 1;


noiseLevel = 5;


dirPath = 'dataset/DeformedMesh';

DataCell = readDir (dirPath);

if DIM == 2
    for i = 1 : size(DataCell, 2)
        A = DataCell{i};
        DataCell{i} = A([1,2], :);
    end
end


InputData = DataCell{1};


% original formulation
ControlPoints = findControlPoints (InputData, dimFeatureSpace, DIM);

IntrisicMatr = TPSWarpModel (ControlPoints, InternalSmoothParam, DIM);

BendingEnergyMatrix = 8 * pi * IntrisicMatr (1:size(IntrisicMatr, 2), :)

fprintf (2, '\nEigen values of the bending enery matrix \n');

eigBendingEnergy = eig ((BendingEnergyMatrix + BendingEnergyMatrix')/2)

ImageCloud = ControlPoints;

wa = IntrisicMatr * ImageCloud'

energy_WEW =  ImageCloud * BendingEnergyMatrix * ImageCloud'


% feature space formulation Warp = W^T * v(D)


Kw = zeros (dimFeatureSpace, size(InputData, 2));

for i = 1 : size (InputData, 2)
    
    Kw(:,i) = fMapping (InputData(:,i), ControlPoints, IntrisicMatr, DIM);
    
end

Coefficient = inv( Kw * Kw' + SmoothParam * BendingEnergyMatrix ) * Kw;


% identical target cloud

fprintf (2, '\nCase of Identitcal target cloud\n');

targetCloud = DataCell{1};

Weights = Coefficient * targetCloud';

ControlPoints__Weights_Difference = [ControlPoints', Weights,  sqrt(sum((ControlPoints' - Weights).^2, 2)) ]

EstimateCloud = Weights' * Kw;

plotCloud (InputData, targetCloud, EstimateCloud, 'Identical Target Cloud', DIM);

energy_WEW =  Weights' * BendingEnergyMatrix * Weights


% deformed target cloud
% fprintf (2, '\nCase of Deformed target cloud\n');
% 
% targetCloud = noisyData (InputData, noiseLevel);
% 
% Weights = Coefficient * targetCloud';
% 
% ControlPoints__Weights_Difference = [ControlPoints', Weights, sqrt(sum((ControlPoints' - Weights).^2, 2))]
% 
% EstimateCloud = Weights' * Kw;
% 
% plotCloud (InputData, targetCloud, EstimateCloud, 'Deformed Target Cloud', DIM);
% 
% energy_WEW =  Weights' * BendingEnergyMatrix * Weights












        % Thin-Plate-Spline basis function for 3D case
        % it's a radial distance function
        % dp = the differnece of two points
        function dist = TPS_BasisFunc ( dp, DIM )
            
            if (DIM == 2)
               
                rr = dp(1) * dp(1) + dp(2) * dp(2);
                
                dist = rr * log(rr);
               
            elseif (DIM == 3)
                
                dist = - norm (dp);
               
            end
            
        end


        



function ControlPoints = findControlPoints (Data, dimFeatureSpace, DIM)

                    min_xyz = min (Data, [], 2);
                    max_xyz = max (Data, [], 2);
                    spn_xyz = max_xyz - min_xyz;
                    
                    if (DIM==2)
                        
                        ControlPoints = [ min_xyz(1) + rand(1, dimFeatureSpace) * spn_xyz(1);
                                                    min_xyz(2) + rand(1, dimFeatureSpace) * spn_xyz(2); ];
                        
                    elseif (DIM ==3)
  
                        ControlPoints = [ min_xyz(1) + rand(1, dimFeatureSpace) * spn_xyz(1);
                                                    min_xyz(2) + rand(1, dimFeatureSpace) * spn_xyz(2);
                                                    min_xyz(3) + rand(1, dimFeatureSpace) * spn_xyz(3); ];
                        
                    end

end



        function plotCloud (Cloud1, Cloud2, Cloud3, figTitle, DIM)
            
            fig = figure;
            fig.Position = [50, 600, 1200, 400];
            
            
            subplot(1,3,1)
            
            if (DIM == 2)
                
                plot(Cloud1(1,:), Cloud1(2,:), 'r.');
                axis equal;
                title ('Source');
                
                subplot(1,3,2)
                plot(Cloud2(1,:), Cloud2(2,:), 'b.');
                axis equal;
                title ('Target');
                
                subplot(1,3,3)
                plot(Cloud3(1,:), Cloud3(2,:), 'g.');
                axis equal;
                title ('Estiamtion');
                
            elseif (DIM == 3)
                
                plot3(Cloud1(1,:), Cloud1(2,:), Cloud1(3,:), 'r.');
                axis equal;
                title ('Source');
                
                subplot(1,3,2)
                plot3(Cloud2(1,:), Cloud2(2,:), Cloud2(3,:), 'b.');
                axis equal;
                title ('Target');
                
                subplot(1,3,3)
                plot3(Cloud3(1,:), Cloud3(2,:), Cloud3(3,:), 'g.');
                axis equal;
                title ('Estiamtion');
                
            end

            sgtitle(figTitle);
            
        end
        



       function [IntrisicMatr, K, InvK, C] = TPSWarpModel (ControlPoints, InternalSmoothParam, DIM)
            
                    % set control points of TPS
                    
                    dimFeatureSpace = size (ControlPoints, 2);
  
                    % compute the intrisic matrix of TPS
                    
                    K = zeros (dimFeatureSpace, dimFeatureSpace);
                    
                    for it_r = 1 : dimFeatureSpace
                        
                        for it_c = it_r+1 : dimFeatureSpace
                            
                            K (it_r, it_c) = TPS_BasisFunc( ControlPoints(:, it_r) - ControlPoints(:, it_c), DIM );
                            
                        end
                        
                    end
                    
                    K = K + K';
                    
                    for it_d = 1 : dimFeatureSpace
                        
                        K (it_d, it_d) = InternalSmoothParam;
                        
                    end
                    
                    InvK = inv (K);
                    
                    C = [ ControlPoints; ones(1, dimFeatureSpace) ];
                    
                    H = inv(C * InvK * C') * C * InvK;
                    
                    IntrisicMatr = [ InvK - InvK * C' * H;  H ];
                    
                    
                    fprintf (2, '\nEigen values of  K_lambda, i.e., Kernel matrix \n');
                    
                    eigK = eig(K)

                    fprintf (2, '\nEigen values of  P * K * P^T \n');
                    
                    eigPKP = eig(C * K * C')
                    
                    fprintf (2, '\nEigen values of  P * InvK * P^T \n');
                    
                    eigPInvKP = eig(C * InvK * C')
 
        end
        
        
        
        % mapping R^3 -> R^l : point to feature space
        function fImage = fMapping (pVal, ControlPoints, IntrisicMatr, DIM)
                           
                dimFeatureSpace = size (ControlPoints, 2);
                
                nlv = ones (dimFeatureSpace + DIM + 1, 1);
                
                dmatr = ControlPoints - pVal * ones (1, dimFeatureSpace);
                
                for k = 1 : dimFeatureSpace
                    
                    nlv (k) = TPS_BasisFunc (dmatr(:, k), DIM);
                    
                end
                
                nlv (dimFeatureSpace+1 : dimFeatureSpace+DIM) = pVal;
                
                fImage =IntrisicMatr' * nlv;
                
        end
        
        

            

        

        
        
          
        function DataCell = readDir (dirPath)
            
            fileHandles = dir ([dirPath, '/*.', 'txt']);
            
            dataMatr = [];
                        
            for ii = 1 : size (fileHandles, 1)
                
                fid = fopen ([dirPath, '/', fileHandles(ii).name], 'r');
                
                lstr = fgetl (fid);
                
                pcnt = 0;
                
                while ischar (lstr)
                    
                    if (size (lstr, 2) > 0)
                        
                        pcnt = pcnt + 1;
                        
                        dataMatr(:, pcnt) = str2double(split(lstr));
                        
                    end
                    
                    lstr = fgetl (fid);
                    
                end
                
                DataCell{ii} = dataMatr;
                
                fclose (fid);
                
            end
            
        end
        
        
        
        function DataNoisy = noisyData (Data, noiselevel)
            
            DataNoisy = Data;
        
            dim = size(Data, 1);
            
            numPoints = size(Data, 2);
            
            for i = 1 : numPoints
                
                DataNoisy(:, i) = Data(:, i) + noiselevel * randn(dim,1);
                
            end

        end
        
        
        
            
        