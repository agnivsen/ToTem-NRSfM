clear all;
close all;
clc;

%-%
addpath('../DefGPA', '../');

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end

flagJointlyOptimizeLamdaRigidTransformation = false;
%-%


% File structure
%
% -- code
% -- dataset/LiverSimu/seq1    
%
filepath = '../../dataset/LiverSimu/seq1/';



% read vertex file
for ii = 1 : 10
    tmp = [filepath, 'deformed_', num2str(ii), '.off'];
    A = readmatrix(tmp, 'FileType','text', 'Range', [3, 1, 2004, 3]);
    DataCell{ii} = (A)';
end

% read visibility of vertex
for ii = 1 : 10
    tmp = [filepath, 'index_', num2str(ii), '.txt'];
    v = readmatrix(tmp, 'FileType','text');
    
    vis = false(1, 2002); vis(v) = true;
    % there will be errors if we use partial observability here
    vis(1:end) = true;  % assume all points are observed.
    vis(10*ii:20*ii) = false;

    VisVecCell{ii} = vis;
end




fprintf(1, " ---- runing GPA_RIGID ----\n");
GPA_RIGID = RigidGPA;
GPA_RIGID.addDataByArray(DataCell, VisVecCell);
GPA_RIGID.run();
GPA_RIGID



fprintf(1, " ---- runing GPA_AFFINE ----\n");
GPA_AFFINE = DefGPA('AFFINE');
GPA_AFFINE.dimFeatureSpace = 0;
GPA_AFFINE.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_AFFINE.addDataByArray(DataCell, VisVecCell);
GPA_AFFINE.run (0, flagJointlyOptimizeLamdaRigidTransformation);
GPA_AFFINE.rsdError();
GPA_AFFINE




% fprintf(1, " ---- runing GPA_TPS ----\n");
% GPA_TPS = DefGPA('TPS');
% GPA_TPS.dimFeatureSpace =5^3;
% GPA_TPS.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
% GPA_TPS.addDataByArray(DataCell, VisVecCell);
% GPA_TPS.run (0.01, flagJointlyOptimizeLamdaRigidTransformation);
% GPA_TPS.rsdError();
% GPA_TPS




% fprintf(1, " ---- runing GPA_KERNEL ----\n");
% GPA_KERNEL = KernelGPA;
% GPA_KERNEL.smoothParam = 0.1;
% GPA_KERNEL.quantilePercentParam= 0.2;
% GPA_KERNEL.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
% GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
% GPA_KERNEL.run ();
% GPA_KERNEL.rsdError();
% GPA_KERNEL





%%
% transform the point-clouds to the reference frame
% for ii = 1 : length(DataCell)
%     data = GPA_RIGID.transformPoints(DataCell{ii}, ii);
%     transformedDataCell{ii} = data;
% end
% hFig100 = figure(100);
% KernelGPA.VisualizePointCloudsVariation (transformedDataCell, VisVecCell)
% title('Rigid'); 
% axis tight equal;
% pause(0.2);



% transform the point-clouds to the reference frame
% for ii = 1 : length(DataCell)
%     data = GPA_TPS.transformPoints(DataCell{ii}, ii);
%     transformedDataCell{ii} = data;
% end
% hFig200 = figure(200);
% KernelGPA.VisualizePointCloudsVariation (transformedDataCell, VisVecCell)
% title('TPS'); 
% axis tight equal;
% pause(0.2);



% transform the point-clouds to the reference frame
% for ii = 1 : length(DataCell)
%     data = GPA_KERNEL.transformPoints(DataCell{ii}, ii);
%     transformedDataCell{ii} = data;
% end
% hFig300 = figure(300);
% KernelGPA.VisualizePointCloudsVariation (transformedDataCell, VisVecCell)
% title('Kernel'); 
% axis tight equal;
% pause(0.2);


