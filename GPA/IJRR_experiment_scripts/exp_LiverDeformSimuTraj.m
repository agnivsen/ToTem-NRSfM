clear all;
close all;
clc;



%-%
addpath('../DefGPA', '../');

fdir = '../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end

flagJointlyOptimizeLamdaRigidTransformation = true;
%-%


% File structure
%
% -- code
% -- dataset/LiverSimu/seq1    
%
filepath = '../../dataset/LiverSimu/downsampled_noise_5mm/';


% idx = 1:4:2002;

% read vertex file
for ii = 0 : 59
    tmp = [filepath, 'full_visible_', num2str(ii), '.pcd'];
    cloud = pcread(tmp);
    A = cloud.Location;
    % A = readmatrix(tmp, 'FileType','text', 'Range', [3, 1, 2004, 3]);
    DataCell{ii+1} = (A)';
end

% read visibility of vertex
for ii = 1 : 60
    % tmp = [filepath, 'index_', num2str(ii), '.txt'];
    % v = dlmread(tmp);
    % v = readmatrix(tmp, 'FileType','text');
    vis = true(1, size(DataCell{ii}, 2)); 
%     vis(randi(size(DataCell{ii}, 2), 1, 500)) = false;
    % vis(idx) = true;
    % there will be errors if we use partial observability here
%     vis(1:end) = true;  % assume all points are observed.

    VisVecCell{ii} = vis;
end




fprintf(1, ' ---- runing GPA_RIGID ----\n');
GPA_RIGID = RigidGPA;
GPA_RIGID.addDataByArray(DataCell, VisVecCell);
GPA_RIGID.run();
GPA_RIGID

%save_traj(GPA_RIGID.mPoses, '../results/poses/gpa_rigid_poses.txt');


fprintf(1, ' ---- runing GPA_AFFINE ----\n');
GPA_AFFINE = DefGPA('AFFINE');
GPA_AFFINE.dimFeatureSpace = 0;
GPA_AFFINE.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_AFFINE.addDataByArray(DataCell, VisVecCell);
GPA_AFFINE.run (0, flagJointlyOptimizeLamdaRigidTransformation);
GPA_AFFINE.rsdError();
GPA_AFFINE
%save_traj(GPA_AFFINE.mPoses, '../results/poses/gpa_affine_poses.txt');




fprintf(1, ' ---- runing GPA_TPS ----\n');
GPA_TPS = DefGPA('TPS');
GPA_TPS.dimFeatureSpace =5^3;
GPA_TPS.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_TPS.addDataByArray(DataCell, VisVecCell);
GPA_TPS.run (0.01, flagJointlyOptimizeLamdaRigidTransformation);
GPA_TPS.rsdError();
GPA_TPS
%save_traj(GPA_TPS.mPoses, '../results/poses/gpa_tps_poses.txt');




fprintf(1, ' ---- runing GPA_KERNEL ----\n');
GPA_KERNEL = KernelGPA;
GPA_KERNEL.smoothParam = 1;
GPA_KERNEL.quantilePercentParam= 0.25;
GPA_KERNEL.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
GPA_KERNEL.run ();
GPA_KERNEL.rsdError();
GPA_KERNEL
%save_traj(GPA_TPS.mPoses, '../results/poses/gpa_kernel_poses.txt');




%%
% transform the point-clouds to the reference frame
for ii = 1 : length(DataCell)
    data = GPA_RIGID.transformPoints(DataCell{ii}, ii);
    transformedDataCell{ii} = data;
    cloud = pointCloud(data');
    %pcwrite(cloud, '../results/deformed_clouds/gpa_rigid_'+string(ii)+'.ply')
end
hFig100 = figure(100);
KernelGPA.VisualizePointCloudsVariation (transformedDataCell, VisVecCell)
title('Rigid'); 
axis tight equal;
pause(0.2);



% transform the point-clouds to the reference frame
for ii = 1 : length(DataCell)
    data = GPA_TPS.transformPoints(DataCell{ii}, ii);
    transformedDataCell{ii} = data;
    cloud = pointCloud(data');
    %pcwrite(cloud, '../results/deformed_clouds/gpa_tps_'+string(ii)+'.ply')
end
hFig200 = figure(200);
KernelGPA.VisualizePointCloudsVariation (transformedDataCell, VisVecCell)
title('TPS'); 
axis tight equal;
pause(0.2);



% transform the point-clouds to the reference frame
for ii = 1 : length(DataCell)
    data = GPA_KERNEL.transformPoints(DataCell{ii}, ii);
    transformedDataCell{ii} = data;
    cloud = pointCloud(data');
    %pcwrite(cloud, '../results/deformed_clouds/gpa_kernel_'+string(ii)+'.ply')
end
hFig300 = figure(300);
KernelGPA.VisualizePointCloudsVariation (transformedDataCell, VisVecCell)
title('Kernel'); 
axis tight equal;
pause(0.2);


% transform the point-clouds to the reference frame
for ii = 1 : length(DataCell)
    data = GPA_AFFINE.transformPoints(DataCell{ii}, ii);
    transformedDataCell{ii} = data;
    cloud = pointCloud(data');
    %pcwrite(cloud, '../results/deformed_clouds/gpa_affine_'+string(ii)+'.ply')
end
hFig400 = figure(400);
KernelGPA.VisualizePointCloudsVariation (transformedDataCell, VisVecCell)
title('Affine'); 
axis tight equal;
pause(0.2);


function [] = save_traj(traj, path)
    aligned_traj = {};
    for i = 1:size(traj,2)
        pi = [traj{i}; 0, 0, 0, 1];
        aligned_traj{i} = pi;
    end
    aligned_traj3x4 = {};
    for i = 1:size(traj,2)
        aligned_traj3x4{i} = inv(aligned_traj{1})*aligned_traj{i};
        aligned_traj3x4{i}(4,:) = [];
    end
    poses = [];
    for i= 1:size(traj,2)
        pi = reshape(aligned_traj3x4{i}, 1, 12);
        poses = [poses; pi];
    end
    fid=fopen(path,'w');
    for i = 1:size(traj,2)
        fprintf(fid,'%.3f,',poses(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end