clear all
close all
clc


addpath('../', '../DefGPA');

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end


figformat = 'pdf';


dim  =3;

TPS_SmoothParam = 0.1;

MC_Trials = 20;



%SigmaArray = [10:10:150];

%SigmaArray = [10, 50, 100, 150, 200];

SigmaArray = 0 : 0.2 : 1;

visualize_sigma = 1;


colors = parula(6);



RMSE_rot_cell = cell(1, length(SigmaArray));
RMSE_pos_cell = cell(1, length(SigmaArray));



for noise_cnt = 1 : length(SigmaArray)


RMSE_rot = [];
RMSE_pos = [];

deform_sigma =  SigmaArray(noise_cnt);

for mc_cnt = 1 : MC_Trials


visualize_example = true;
if (visualize_example)
    deform_sigma = visualize_sigma;
end


simulator = DeformingPointCloudSimulator(80 , 100);
simulator.sigma_deform = deform_sigma;

[DataCell, VisVecCell] =  simulator.SimulateMeasurements ();


vis_all = zeros(1, length(VisVecCell{1}));
for ii = 1 : length(VisVecCell)
    vis_all = vis_all + VisVecCell{ii};
end
vis_all



GT_Landmarks = simulator.GT_InitialLandmarks;
GT_Poses = simulator.GT_Poses;

if (0)
figure("Name", "Ground_truth", "Position", [0,0,600, 500]);
[h0, hl0] = VisualizePoses (GT_Poses, GT_Landmarks, eye(4));  h0.Color = 'k';  hl0.MarkerEdgeColor = 'k';
control_points = simulator.TPS.ControlPoints;
hold on
scatter3(control_points(1, :), control_points(2, :), control_points(3, :), 'bo', 'filled');
hold off;
axis equal; view(3);

figure("Name", "PointCloudSample Partial", "Position", [0,800,600, 500]); hold on;
[h0, hl0] =VisualizePoses (simulator.GT_Poses, simulator.GT_DeformedLandmarks{1}, eye(4));  h0.Color = 'r';  hl0.MarkerEdgeColor = 'r';
%[h0, hl0] =VisualizePoses (simulator.GT_Poses, simulator.GT_DeformedLandmarks{2}, eye(4));  h0.Color = 'g';  hl0.MarkerEdgeColor = 'g';
%[h0, hl0] =VisualizePoses (simulator.GT_Poses, simulator.GT_DeformedLandmarks{3}, eye(4));  h0.Color = 'b';  hl0.MarkerEdgeColor = 'b';
hold off; axis equal; view(3);

figure("Name", "PointCloudSample All", "Position", [800,800,600, 500]);
for ii = 1 : length(DataCell)
    [h0, hl0] =VisualizePoses (simulator.GT_Poses, simulator.GT_DeformedLandmarks{ii}, eye(4));  h0.Color = 'b';  hl0.MarkerEdgeColor = 'b';
end
axis equal; view(3);

return;
end


flagJointlyOptimizeLamdaRigidTransformation = true;


GPA_RIGID = RigidGPA;
GPA_RIGID.addDataByArray(DataCell, VisVecCell);
GPA_RIGID.flagEliminateReferenceShapeFirst = false;
GPA_RIGID.run();
GPA_RIGID



GPA_AFFINE = DefGPA('AFFINE');
GPA_AFFINE.dimFeatureSpace = 0;
GPA_AFFINE.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_AFFINE.addDataByArray(DataCell, VisVecCell);
GPA_AFFINE.run (0, flagJointlyOptimizeLamdaRigidTransformation);
GPA_AFFINE.rsdError();
GPA_AFFINE



GPA_TPS5 = DefGPA('TPS');
GPA_TPS5.dimFeatureSpace =5^dim;
GPA_TPS5.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_TPS5.addDataByArray(DataCell, VisVecCell);
GPA_TPS5.run (TPS_SmoothParam, flagJointlyOptimizeLamdaRigidTransformation);
GPA_TPS5.rsdError();
GPA_TPS5



GPA_TPS7 = DefGPA('TPS');
GPA_TPS7.dimFeatureSpace = 7^dim;
GPA_TPS7.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_TPS7.addDataByArray(DataCell, VisVecCell);
GPA_TPS7.run (TPS_SmoothParam, flagJointlyOptimizeLamdaRigidTransformation);
GPA_TPS7.rsdError();
GPA_TPS7



GPA_KERNEL = KernelGPA;
GPA_KERNEL.smoothParam = 0.05;
GPA_KERNEL.quantilePercentParam= 0.2;
GPA_KERNEL.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
GPA_KERNEL.run ();
GPA_KERNEL.rsdError();
GPA_KERNEL



[R, t] = EuclideanProcrustes (GPA_RIGID.mShape, GT_Landmarks);   T1 =[R, t];
[R, t] = EuclideanProcrustes (GPA_AFFINE.mShape, GT_Landmarks);  T2 = [R, t];
[R, t] = EuclideanProcrustes (GPA_TPS5.mShape, GT_Landmarks);  T3 = [R, t];
[R, t] = EuclideanProcrustes (GPA_TPS7.mShape, GT_Landmarks);  T4 = [R, t];
[R, t] = EuclideanProcrustes (GPA_KERNEL.mShape, GT_Landmarks);  T5 = [R, t];



[RMSE_rot_1, RMSE_pos_1] = TrajectoryError.RelativePoseError (GPA_RIGID.mPoses, GT_Poses);
[RMSE_rot_2, RMSE_pos_2] = TrajectoryError.RelativePoseError (GPA_AFFINE.mPoses, GT_Poses);
[RMSE_rot_3, RMSE_pos_3] = TrajectoryError.RelativePoseError (GPA_TPS5.mPoses, GT_Poses);
[RMSE_rot_4, RMSE_pos_4] = TrajectoryError.RelativePoseError (GPA_TPS7.mPoses, GT_Poses);
[RMSE_rot_5, RMSE_pos_5] = TrajectoryError.RelativePoseError (GPA_KERNEL.mPoses, GT_Poses);


tmp_rot = [RMSE_rot_1, RMSE_rot_2, RMSE_rot_3, RMSE_rot_4, RMSE_rot_5]
tmp_pos = [RMSE_pos_1, RMSE_pos_2, RMSE_pos_3, RMSE_pos_4, RMSE_pos_5]

RMSE_rot = [RMSE_rot;  tmp_rot];
RMSE_pos = [RMSE_pos;  tmp_pos];

if (visualize_example && (mc_cnt == 1) &&(deform_sigma == visualize_sigma))
    
    hFig10 = figure(10); 
    tfig10 = tiledlayout(1, 1);
    ax = nexttile;
    hold on;
    [h0, hl0] = VisualizePoses (GT_Poses, GT_Landmarks, eye(4));  h0.Color = 'r';  hl0.MarkerEdgeColor = 'r';
    [h1, hl1] = VisualizePoses (GPA_RIGID.mPoses, GPA_RIGID.mShape, T1);  h1.Color = colors(1, :); hl1.MarkerEdgeColor = colors(1, :); hl1.MarkerFaceColor = colors(1, :);
    [h2, hl2] = VisualizePoses (GPA_AFFINE.mPoses, GPA_AFFINE.mShape, T2); h2.Color = colors(2, :);  hl2.MarkerEdgeColor = colors(2, :);  hl2.MarkerFaceColor = colors(2, :);
    [h3, hl3] = VisualizePoses (GPA_TPS5.mPoses, GPA_TPS5.mShape,  T3); h3.Color = colors(3, :);  hl3.MarkerEdgeColor = colors(3, :); hl3.MarkerFaceColor = colors(3, :);
    [h4, hl4] = VisualizePoses (GPA_TPS7.mPoses, GPA_TPS7.mShape,  T4); h4.Color = colors(4, :);  hl4.MarkerEdgeColor = colors(4, :); hl4.MarkerFaceColor = colors(4, :);
    [h5, hl5] = VisualizePoses (GPA_KERNEL.mPoses, GPA_KERNEL.mShape, T5); h5.Color = colors(5, :);  hl5.MarkerEdgeColor = colors(5, :); hl5.MarkerFaceColor = colors(5, :);
    axis equal tight; axis on;    box off;
    ax.FontSize = 8;
    ax.LineWidth = 0.2;
    ax.YAxis.Visible = 'off'; 
    ax.ZAxis.Visible = 'off';
    ax.XAxis.TickDirection = "in";
    ax.XAxis.TickLength = [0.01, 0.01];
    view(0, 10);

    lgd = legend([hl0, hl1, hl2, hl3, hl4, hl5],  {'Ground-Truth',  'EUC', 'AFF', 'TPS(5)', 'TPS(7)', 'Kernel'},  ...
         'FontSize', 9, 'Orientation', 'vertical',  'box', 'off');
    lgd.Layout.Tile = 'east';

    hFig10.Position = [0 , 600,  550, 250];
    exportgraphics(hFig10, [fdir, 'PoseEstimationError_ExampleVisualization', '.pdf'], 'ContentType','vector');

    fig100 = figure("Name", "Fitness of each methods", "Position", [0, 100, 900, 1200]);
    tfig100 = tiledlayout(3, 2);
    nexttile; VisualizeRegistration (GPA_RIGID);  view(3)
    nexttile; VisualizeRegistration (GPA_AFFINE);  view(3)
    nexttile; VisualizeRegistration (GPA_TPS5);  view(3)
    nexttile; VisualizeRegistration (GPA_TPS7);  view(3)
    nexttile; VisualizeRegistration (GPA_KERNEL);  view(3)
    nexttile; VisualizePoses (GT_Poses, GT_Landmarks, eye(4));  view(3)
    return;
    
end


end

RMSE_rot_cell{noise_cnt} = RMSE_rot;
RMSE_pos_cell{noise_cnt} = RMSE_pos;



end


return;

xvalues = SigmaArray;

[rot_rigid, rot_affine, rot_tps5, rot_tps7, rot_kernel] = ProcessDataCell(RMSE_rot_cell);
[pos_rigid, pos_affine, pos_tps5, pos_tps7, pos_kernel] = ProcessDataCell(RMSE_pos_cell);

save('PoseEstimationErrorMonteCarlo', 'RMSE_rot_cell', 'RMSE_pos_cell', 'SigmaArray', 'xvalues', ...
    "rot_rigid", 'rot_affine', "rot_tps5", 'rot_tps7', 'rot_kernel', ...
    "pos_rigid", "pos_affine", "pos_tps5", "pos_tps7", "pos_kernel");



%%% -------------- visualization -------------- %%%

labelfontSize = 8;

hFig1 = figure(1);
hFig1.Position = [0 , 600,  500, 250];

tfig1 = tiledlayout(1, 2, 'TileSpacing','compact');


ax2 = nexttile; 
ax2.TickDir = 'in';
ax2.FontSize = 7;
hold on;
pos_RMSE_mean = [mean(pos_rigid); mean(pos_affine); mean(pos_tps5); mean(pos_tps7); mean(pos_kernel)];
hb2 = bar(xvalues, pos_RMSE_mean', 'FaceColor','flat',  'EdgeColor', 'none');
for i = 1 : 5
    xx(:, i) = hb2(1, i).XEndPoints';
    hb2(i).FaceColor = colors(i, :);
end
pos_RMSE_std = [std(pos_rigid); std(pos_affine); std(pos_tps5); std(pos_tps7); std(pos_kernel)];
errorbar(xx, pos_RMSE_mean', pos_RMSE_std', 'LineStyle', 'none', 'color', 'c', 'LineWidth', 1.0);
hold off; box on; grid on; axis tight;
xlabel('$\sigma_{\mathrm{deform}}$', 'Interpreter','latex', 'FontSize', labelfontSize+2);
ylabel('translational RPE', 'FontSize', labelfontSize);
xticks(xvalues);


ax1 = nexttile;
ax1.TickDir = 'in';
ax1.FontSize = 7;
hold on;
rot_RMSE_mean = [mean(rot_rigid); mean(rot_affine); mean(rot_tps5); mean(rot_tps7); mean(rot_kernel)] ;
hb1 = bar(xvalues, rot_RMSE_mean', 'FaceColor','flat', 'EdgeColor', 'none');
rot_RMSE_std = [std(rot_rigid); std(rot_affine); std(rot_tps5); std(rot_tps7); std(rot_kernel)] ;
for i = 1 : 5
    xx(:, i) = hb1(1, i).XEndPoints';
    hb1(i).FaceColor = colors(i, :);
end
errorbar(xx, rot_RMSE_mean', rot_RMSE_std', 'LineStyle', 'none', 'color', 'c', 'LineWidth', 1.0);
hold off; box on; grid on; axis tight;
xlabel('$\sigma_{\mathrm{deform}}$', 'Interpreter','latex', 'FontSize', labelfontSize+2);
ylabel('rotational RPE', 'FontSize', labelfontSize);
xticks(xvalues);

lgd = legend(hb1,  {'EUC', 'AFF', 'TPS(5)', 'TPS(7)', 'Kernel'},  'FontSize', 7, 'Orientation', 'Horizontal', 'box', 'off');
lgd.Layout.Tile = 'north';


exportgraphics(tfig1, [fdir, 'PoseEstimationError_MC', '.pdf'], 'ContentType','vector');






function [arr1, arr2, arr3, arr4, arr5] = ProcessDataCell(mcStatisticsCell)
    
    arr1 = [];
    arr2 = [];
    arr3 = [];
    arr4 = [];
    arr5 = [];

    for ii = 1 : length(mcStatisticsCell)
        vals = mcStatisticsCell{ii};
        arr1 = [arr1,  vals(:, 1)];
        arr2 = [arr2,  vals(:, 2)];
        arr3 = [arr3,  vals(:, 3)];
        arr4 = [arr4,  vals(:, 4)];
        arr5 = [arr5,  vals(:, 5)];
    end
end






function VisualizeRegistration (GPA_Handle) 
    rshape = GPA_Handle.mShape;
    num_poses = GPA_Handle.numPointClouds;
    PoseArray = GPA_Handle.mPoses;

    hold on;
    scatter3 ( rshape(1, :), rshape(2, :), rshape(3, :), 'ro');
    for ii = 1 : num_poses
        pts = GPA_Handle.PointClouds(ii).Data;
        vis = GPA_Handle.PointClouds(ii).Vis;
        pts_vis = pts(:, logical(vis));

        transform_pts = GPA_Handle.transformPoints(pts_vis, ii);        
        scatter3( transform_pts(1, :), transform_pts(2, :), transform_pts(3, :), '*');
    end
    VisualizePoses (PoseArray, rshape);
    hold off;

end




function [hp, hf] = VisualizePoses (PoseArray,  landmarks,  GaugeT)

if ~exist('GaugeT', 'var')
    GaugeT = eye(4);
end

nsize = length(PoseArray);

pose = PoseArray{1};
if (size(pose, 1) == size(pose, 2))
    ndim = size(pose, 1) - 1;
end
if (size(pose, 1) + 1 == size(pose, 2))
    ndim = size(pose, 1);
end

GaugeRot = GaugeT(1:ndim, 1:ndim);
GaugePos = GaugeT(1:ndim, 1+ndim);

hold on;

trajectory_position = [];

for ii = 1 : nsize
    pose = PoseArray{ii};
    R = pose(1:ndim, 1:ndim);
    t = pose(1:ndim, 1+ndim);
    
    roation = GaugeRot*R;
    position = GaugeRot*t + GaugePos;
    
    trajectory_position = [trajectory_position,  position];
end


if ndim == 2
    hp = plot(trajectory_position(1, :), trajectory_position(2, :), 'r.-');
end

if ndim == 3
    hp = plot3(trajectory_position(1, :), trajectory_position(2, :), trajectory_position(3, :), 'r.-');
    view(3);
end


if exist('landmarks', 'var') && ~isempty(landmarks)
    landmarks = GaugeRot * landmarks + GaugePos;
    if size(landmarks, 1) == 2
        hf = scatter(landmarks(1, :), landmarks(2, :), 4, 'o', 'filled');
    end    
    if size(landmarks, 1) == 3
        hf = scatter3(landmarks(1, :), landmarks(2, :), landmarks(3, :), 4, 'o', 'filled');
        view(3);
    end
end

hold off;    

end






% R * D1 + t = D2
function [R, t] = EuclideanProcrustes (D1, D2)

meanD1 = mean(D1, 2);

meanD2 = mean(D2, 2);

M = (D1 - meanD1) * (D2 - meanD2)';

[U, ~, V] = svd(M);

R = V * U';

if det(R) < 0
    
    fprintf(2, '\n There exists reflection between solutions! \n');
    
end

t =  meanD2 - R * meanD1;

end



 
 