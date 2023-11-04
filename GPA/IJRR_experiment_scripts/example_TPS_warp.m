clear all;
clc;


addpath('../DefGPA', '../');


% input_point_cloud,     3 * N matrix
% deformed_point_cloud   3 * N matrix
deformed_point_cloud = UseTPS_Method1(input_point_cloud, controlPointsInitialPosition, controlPointsTargetPosition);




% defomation is controlld by control points,
% controlPointsInitialPosition   ---->  controlPointsTargetPosition
function deformed_point_cloud = UseTPS_Method1(input_point_cloud, controlPointsInitialPosition, controlPointsTargetPosition)
% input_point_cloud,     3 * N matrix
% deformed_point_cloud   3 * N matrix

deformed_point_cloud = zeros(size(input_point_cloud, 1), size(input_point_cloud, 2));

InternalSmoothParam = 0; % diagonal elements of TPS kernel matrix. default = 0
IntrisicMatr = ThinPlateSpline.setTPSIntrisicMatr (controlPointsInitialPosition, InternalSmoothParam);

for ii = 1 : size(input_point_cloud, 2)
  nlv = ThinPlateSpline.getEtaTPS (controlPointsInitialPosition, input_point_cloud(:, ii));
  deformed_point_cloud(:, ii) = controlPointsTargetPosition * IntrisicMatr' * nlv;
end

end


% deformation is controlled by a single scalar sigma_deform
% control points are chosen automatically, 5 points along each principal axis.
function deformed_point_cloud = UseTPS_Method2(input_point_cloud, sigma_deform)
% input_point_cloud,     3 * N matrix
% deformed_point_cloud   3 * N matrix

deformed_point_cloud = zeros(size(input_point_cloud, 1), size(input_point_cloud, 2));

dimFeatureSpace = 5^size(input_point_cloud, 1);   % 5 points along each principal axis
controlPointsInitialPosition =  ThinPlateSpline.setControlPoints (input_point_cloud, dimFeatureSpace);

InternalSmoothParam = 0; % diagonal elements of TPS kernel matrix. default = 0
IntrisicMatr = ThinPlateSpline.setTPSIntrisicMatr (controlPointsInitialPosition, InternalSmoothParam);

controlPointsTargetPosition = controlPointsInitialPosition + sigma_deform * max(min(randn(size(controlPointsInitialPosition, 1),  size(controlPointsInitialPosition, 2)), 3), -3);
for ii = 1 : size(input_point_cloud, 2)
  nlv = ThinPlateSpline.getEtaTPS (controlPointsInitialPosition, input_point_cloud(:, ii));
  deformed_point_cloud(:, ii) = controlPointsTargetPosition * IntrisicMatr' * nlv;
end

end