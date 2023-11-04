clear all
close all
clc

addpath('DefGPA');


% choose one dataset here
dchoice = 7;

% ----- Dataset Information ---

datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug', 'TOPACS'};
datasetName = datasetCell{dchoice};
[DataCell, VisVecCell] = readDataset (datasetName);



bestParams = {10, 100, 100, 0.1, 0.01, 0.01,0.01};
dimensions = {2, 2, 2, 3, 3, 3, 3}; 


% ----- Example code on how to use the KernelGPA class ---

KGPA_bestQuantilePercentParams = [0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2];
KGPA_bestSmoothParams = [0.05, 0.05, 0.05, 0.05,  0.05, 0.05, 0.05];


flagJointlyOptimizeLamdaRigidTransformation = true;

% Nfold = inf;
% By setting Nfold to infinity, the method
% run_N_Fold_CrossValidation(Nfold)
% generates Leave-1-out cross-validation.
Nfold = 20;   




GPA_RIGID = RigidGPA;
GPA_RIGID.addDataByArray(DataCell, VisVecCell);
GPA_RIGID.run();
GPA_RIGID
% 
% 
% GPA_AFFINE = DefGPA('AFFINE');
% GPA_AFFINE.dimFeatureSpace = 0;
% GPA_AFFINE.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
% GPA_AFFINE.addDataByArray(DataCell, VisVecCell);
% GPA_AFFINE.run (bestParams{dchoice}, flagJointlyOptimizeLamdaRigidTransformation);
% GPA_AFFINE.rsdError();
% GPA_AFFINE


% GPA_TPS = DefGPA('TPS');
% GPA_TPS.dimFeatureSpace =5^dimensions{dchoice};
% GPA_TPS.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
% GPA_TPS.addDataByArray(DataCell, VisVecCell);
% GPA_TPS.run (bestParams{dchoice}, flagJointlyOptimizeLamdaRigidTransformation);
% GPA_TPS.rsdError();
% GPA_TPS



% GPA_KERNEL = KernelGPA;
% GPA_KERNEL.smoothParam = KGPA_bestSmoothParams(dchoice);
% GPA_KERNEL.quantilePercentParam= KGPA_bestQuantilePercentParams(dchoice);
% GPA_KERNEL.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
% GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
% GPA_KERNEL.run ();
% GPA_KERNEL.rsdError();
% GPA_KERNEL




%%
% transform the point-clouds to the reference frame
for ii = 1 : length(DataCell)
    data = GPA_RIGID.transformPoints(DataCell{ii}, ii);
    transformedDataCell{ii} = data;
end
hFig100 = figure(100);
KernelGPA.VisualizePointCloudsVariation (transformedDataCell, VisVecCell)
title('Rigid'); 
pause(0.2);
% 
% 
% 
% % transform the point-clouds to the reference frame
% for ii = 1 : length(DataCell)
%     data = GPA_AFFINE.transformPoints(DataCell{ii}, ii);
%     transformedDataCell{ii} = data;
% end
% hFig200 = figure(200);
% KernelGPA.VisualizePointCloudsVariation (transformedDataCell, VisVecCell)
% title('AFFINE'); 
% pause(0.2);



% transform the point-clouds to the reference frame
% for ii = 1 : length(DataCell)
%     data = GPA_TPS.transformPoints(DataCell{ii}, ii);
%     transformedDataCell{ii} = data;
% end
% hFig300 = figure(300);
% KernelGPA.VisualizePointCloudsVariation (transformedDataCell, VisVecCell)
% title('TPS'); 
% pause(0.2);



% transform the point-clouds to the reference frame
% for ii = 1 : length(DataCell)
%     data = GPA_KERNEL.transformPoints(DataCell{ii}, ii);
%     transformedDataCell{ii} = data;
% end
% hFig400 = figure(400);
% KernelGPA.VisualizePointCloudsVariation (transformedDataCell, VisVecCell)
% title('Kernel'); 
% pause(0.2);
