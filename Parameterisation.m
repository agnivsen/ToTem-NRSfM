%% Parameterisation
% This code replicates our results from figure 6 of our ToTem  NRSfM 
% paper [1] (pre-print provided in the 'Docs' folder). Note that in this script 
% we only provide qualitative comparison with Ball-Pivot (BP) method using
% the code from [3] which is based on [2]. For Poisson Surface
% Reconstruction (PSR) [4], we used Matlab's inbuilt function
% 'pc2surfacemesh' from [5] which, unfortunuately, is available only with
% Matlab version 2022b and later. The code from [3] is already included in
% the 'Dependencies' folder. We also need a ray-casting method for visualising 
% the reconstructed meshes; we used the method from [6] and has been provided 
% in the 'Dependencies' folder as well.
% -----------------------------------------------------------------------------------
% [1]: Sengupta, Agniva and Adrien Bartoli. "ToTem NRSfM: Object-wise Non-Rigid 
% Structure-from-Motion with a Topological Template" International journal of computer 
% vision, accepted Sep, 2023.
%
% [2]: Bernardini, F., Mittleman, J., Rushmeier, H., Silva, C., & Taubin, G. (1999). The 
% ball-pivoting algorithm for surface reconstruction. IEEE transactions on visualization 
% and computer graphics, 5(4), 349-359.
%
% [3]: Luigi Giaccari (2023). Surface Reconstruction from scattered points cloud (open surfaces) 
% (mathworks.com/matlabcentral/fileexchange/63731-surface-reconstruction-from-scattered-points-cloud-open-surfaces), 
% MATLAB Central File Exchange. Retrieved October 20, 2023.
%
% [4]: Kazhdan, M., Bolitho, M., & Hoppe, H. (2006, June). Poisson surface reconstruction. 
% In Proceedings of the fourth Eurographics symposium on Geometry processing (Vol. 7, p. 0).
%
% [5]: https://fr.mathworks.com/help/lidar/ref/pc2surfacemesh.html
%
% [6]: Vipin Vijayan (2023). Ray casting for deformable triangular 3D meshes 
% (mathworks.com/matlabcentral/fileexchange/41504-ray-casting-for-deformable-triangular-3d-meshes), 
% MATLAB Central File Exchange. Retrieved October 20, 2023.


clear vars; clear all; close all; clear path;

% adding all necessary paths
addpath('./Utils/');
addpath('./Parameterisation/');
addpath('./Dependencies/');
addpath('./Dependencies/gloptipoly3/');
addpath('./Dependencies/SeDuMi_1_3/');
addpath('./Dependencies/CrustOpen/');

% The options for algebraic shapes
TOPOLOGY_PLANE_ = 0;
TOPOLOGY_CYLINDER_ = 1;
TOPOLOGY_SPHERE_ = 2;

% Select the shape you want to test for: 
TOPOLOGY = TOPOLOGY_CYLINDER_; % <=== Change as needed
% Replace with TOPOLOGY_PLANE_, TOPOLOGY_CYLINDER_, or TOPOLOGY_SPHERE_

% Setting simulation and reconstruction options:
 Npts = 100; % Number of points to simulate
 dataWeight = 0.85;% (1 - \kappa_0 - \kappa_1) in equation (33) of [1]
 normalWeight = 0.1499; % \kappa_0 in equation (33) of [1]
 interpolationDensity = 7000; % Resolution of final dense interpolation of algebraic shape (for visualisation)
 maxIteration = 20; % Max. iteration
 noiseMultiplier = 0.15; % Multiplicative noise to the added to simulated pointcloud

 %% Simulation setup:
 if(TOPOLOGY == TOPOLOGY_PLANE_)
     % Generates a random plane
     [pointcloud, gtPose, gtCloud, testCloud, testCloudGT] =... 
     generateRandomPlane(Npts, noiseMultiplier);
 elseif(TOPOLOGY == TOPOLOGY_CYLINDER_)
     % Generates a random cylinder
     [pointcloud, gtPose, gtCloud, testCloud, testCloudGT] =... 
     generateRandomCylinder(Npts, noiseMultiplier);
 elseif(TOPOLOGY == TOPOLOGY_SPHERE_)
     % Generates a random sphere
     [pointcloud, gtPose, gtCloud, testCloud, testCloudGT] =... 
     generateRandomSphere(Npts, noiseMultiplier);
 end
     

 % Randomly picking indices for test-train split
 indices = randperm(Npts);
 indices = indices(1:(Npts/2));
 pointcloudT = pointcloud(indices,:);
 pointcloud(indices,:) = [];

 %% Running PARAMETERISATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For re-using our code, please just modify the following lines with your
%%% own data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
visibility = ones(1, Npts); 
ignorePts = randperm(Npts);
visibility(ignorePts(1:5)) = 0; % artificially hiding 5 random points

 % Invoking the method for 'parameterisation' for algebraic shapes, as described in section 6 of [1]   
 [interpolated_data, T, cleaned_reconstruction] = InitialRefinement(pointcloudT, visibility,... 
     maxIteration, TOPOLOGY, dataWeight, normalWeight, interpolationDensity);
            


 %% Clean-up & visualisation

transformedData = T(1:3, 1:3)*interpolated_data(:, 1:3).' + T(1:3, 4);
transformedData = transformedData.';

sparsePoints = cleaned_reconstruction(:, 1:3);
sparsePoints = T(1:3, 1:3)*sparsePoints.' + T(1:3, 4);
sparsePoints = sparsePoints.';

normals = cleaned_reconstruction(:, 4:6);
normals =  T(1:3, 1:3)*normals.';
normals = -normals.';


%% Visualization
f1 = figure(1);
f1.Position = [5,5,400,400];
p1 = plot3(pointcloudT(:,1).',pointcloudT(:,2).',pointcloudT(:,3).','o'); hold on;
p2 = plot3(pointcloud(:,1).',pointcloud(:,2).',pointcloud(:,3).','o'); hold on;
p1.MarkerFaceColor = [0.9608 0.0941 0.4706]; 
p1.MarkerEdgeColor = [0.9608 0.0941 0.4706];
p2.MarkerFaceColor = [0.0706 0.4824 0.7804]; 
p2.MarkerEdgeColor = [0.0706 0.4824 0.7804];
legend({'Training point-cloud', 'Testing point-cloud'}, 'Interpreter','latex','FontSize',9);
axis on; grid on; axis equal;
xlabel('X', 'Interpreter','latex','FontSize',9);
ylabel('Y', 'Interpreter','latex','FontSize',9);
zlabel('Z', 'Interpreter','latex','FontSize',9);

f1 = figure(2);
f1.Position = [50,100,1200,274];
subplot(1,2,1);
titleString = strcat('IsoWT-IP-P');
 title(titleString,'Interpreter','latex','FontSize',9); hold on;
stopDisplayingNewPlot = true;
pSurf = PlotSurfaceReconstruction(transformedData);
pSurf = pSurf.addSparsePointcloud(sparsePoints, normals);
pSurf.startShowingAxes();
pSurf.addLightToScene([2 2 1]);
pSurf.addLightToScene([0 0 -1]);
pSurf = pSurf.addLightToScene([-0.5 -0.5 2]);
pSurf.setNormalLength(0.25);
pSurf.setMaterial(1);
pSurf.plotReconstructionNow(stopDisplayingNewPlot);
pause(0.1); hold off;

subplot(1,2,2);
titleString = strcat('Ball~Pivot~Reconstruction');
 title(titleString,'Interpreter','latex','FontSize',9); hold on;
stopDisplayingNewPlot = true;
pSurf = PlotSurfaceReconstruction(pointcloudT);
pSurf.startShowingAxes();
pSurf.addLightToScene([2 2 1]);
pSurf.addLightToScene([0 0 -1]);
pSurf.addLightToScene([-0.5 -0.5 2]);
pSurf.setNormalLength(0.1);
pSurf.setMaterial(1);
pSurf.plotReconstructionNow(stopDisplayingNewPlot);
pause(0.1); hold off;


%% generateRandomPlane
% Generates a random plane
function [pointcloud, gtPose, gtCloud, testCloud, testCloudGT] =... 
                                                                        generateRandomPlane(Npts, noiseMultiplier)
    nPts = Npts;
    
    X = random_number_within_range(-2,2,nPts);
    Y = random_number_within_range(-2,2,nPts);
    Z = zeros(1,nPts);
    pointcloud = [X.' Y.' Z.'];
    
    gtCloud = pointcloud;
    
    R = rotx(randi(20)) * roty(randi(20)) * rotz(randi(20));
    t = rand(3,1);
    
    pointcloud = R*pointcloud.' + t;
    pointcloud = pointcloud.';
    
    pointcloud = pointcloud + (noiseMultiplier.*rand(nPts,3));
    
    numTestPts = 4;
    testZ = random_number_within_range(-1 ,1 , numTestPts);    
    testCloudGT = [zeros(numTestPts,1) zeros(numTestPts,1) testZ.'];
    
    testCloud = R*testCloudGT.' + t;
    testCloud = testCloud.';
    
    gtPose = [R.' (-R.'*t); 0 0 0 1];
end

%% generateRandomCylinder
% Generates a random cylinder
function [pointcloud, gtPose, gtCloud, testCloud, testCloudGT] =... 
                                                                        generateRandomCylinder(Npts, noiseMultiplier)
    nPts = Npts;
    theta = random_number_within_range(-pi/2,pi/2, nPts);
    h = random_number_within_range(0,rand*30, nPts);

    radius = range(h)*0.1;
    X = radius.*sin(theta);
    Y = h;
    Z = radius.*cos(theta);
    pointcloud = [X.' Y.' Z.'];
    
    gtCloud = pointcloud;
    
    R = rotx(randi(20))*roty(randi(20))*rotz(randi(20));
    t = rand(3,1) + [0;0;2];
    
    pointcloud = R*pointcloud.' + t;
    pointcloud = pointcloud.';
    
    pointcloud = pointcloud + (noiseMultiplier.*rand(nPts,3));
    
    numTestPts = 4;
    testY = random_number_within_range(-1 ,1 , numTestPts);    
    testCloudGT = [zeros(numTestPts,1) testY.' zeros(numTestPts,1)];
    
    testCloud = R*testCloudGT.' + t;
    testCloud = testCloud.';
    
    gtPose = [R.' (-R.'*t); 0 0 0 1];
end

%% generateRandomSphere
% Generates a random sphere
function [pointcloud, gtPose, gtCloud, testCloud, testCloudGT] =... 
                                                                    generateRandomSphere(Npts, noiseMultiplier)
    nPts = Npts;
    
    T1 = random_number_within_range(-pi/2,pi/2,nPts);
    T2 = random_number_within_range(0,(2*pi),nPts);
    Z = zeros(1,nPts);
    pointcloud = [(-cos(T1).*cos(T2)).' (-sin(T1)).' (-cos(T1).*sin(T2)).'];
    
    gtCloud = pointcloud;
    
    R = rotx(randi(90))*roty(randi(90))*rotz(randi(90));
    t = rand(3,1) + [0;0;1];
    
    pointcloud = R*pointcloud.' + t;
    pointcloud = pointcloud.';
    
    pointcloud = pointcloud + (noiseMultiplier*rand(nPts,3)); % adding noise to pointcloud
    
    testCloudGT = [0 0 0];
    testCloud = t.';
    
    gtPose = [R.' (-R.'*t); 0 0 0 1];
end