clear all
close all
clc


fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end


% include the files in DefGPA directory
addpath('../DefGPA', '../');

datasetName = 'TOPACS';

ffdir = '../../dataset';


datafiledir = '../../dataset/TOPACS/landmarks/';


class_str = "$\\mathcal{C}_{%d}$";


[DataCell, VisVecCell] = readDataset (datasetName, ffdir);




% Kernel-GPA
GPA_KERNEL = KernelGPA;
GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
for testQuantilePercentParam = 0.1 : 0.1 : 1;
    K = GPA_KERNEL.getKernelMatrix(1, testQuantilePercentParam);
    threshold = exp(-3*3);
    offdiag_sparsity = (nnz(K>threshold) - length(K)) / (length(K)*length(K));
    fprintf(1, "q = %f,  off-diagonal-sparsity = %f \n", testQuantilePercentParam, offdiag_sparsity);
end



if(1)

fid = fopen ('../../PaperDraft/figure/exp_TOPACS_validation_table.tex', 'w');
format_spec = " %s  &  %s  &  %s  & %.2f  &  %.2f  &  %.2f  &  \\textbf{%.2f} \\\\ \n";
fprintf(fid, "\\begin{tabular}{@{} ccc|cccc @{}} \n");
fprintf(fid, "\\toprule \n");
fprintf(fid, " training  &  test  &   & Rigid-GPA & Affine-GPA  & TPS-GPA  &  Kernel-PGA  \\\\ \n");

fprintf(fid, "\\midrule \n");

trainingset_indicator = [3];    testset_indicator = [4, 5, 6];
[statistics, hfig] = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator);
fprintf(fid, format_spec, "",  "", "min (mm)", ...
    statistics.EUC.min, statistics.AFF.min, statistics.TPS5.min, statistics.Kernel.min);
fprintf(fid, format_spec, sprintf(class_str, trainingset_indicator),  sprintf(class_str, testset_indicator), "max (mm)", ...
    statistics.EUC.max, statistics.AFF.max, statistics.TPS5.max, statistics.Kernel.max);
fprintf(fid, format_spec, "",  "", "mean (mm)", ...
    statistics.EUC.mean, statistics.AFF.mean, statistics.TPS5.mean, statistics.Kernel.mean);
fprintf(fid, format_spec, "",  "", "RMSE (mm)", ...
    statistics.EUC.RMSE_ve, statistics.AFF.RMSE_ve, statistics.TPS5.RMSE_ve, statistics.Kernel.RMSE_ve);
exportgraphics(hfig, [fdir, 'TOPACS_GPA_Visualization_3_456', '.pdf'], 'ContentType','vector');

fprintf(fid, "\\midrule \n");

trainingset_indicator = [4];    testset_indicator = [3, 5, 6];
[statistics, hfig] = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator);
fprintf(fid, format_spec, "",  "", "min (mm)", ...
    statistics.EUC.min, statistics.AFF.min, statistics.TPS5.min, statistics.Kernel.min);
fprintf(fid, format_spec, sprintf(class_str, trainingset_indicator),  sprintf(class_str, testset_indicator), "max (mm)", ...
    statistics.EUC.max, statistics.AFF.max, statistics.TPS5.max, statistics.Kernel.max);
fprintf(fid, format_spec, "",  "", "mean (mm)", ...
    statistics.EUC.mean, statistics.AFF.mean, statistics.TPS5.mean, statistics.Kernel.mean);
fprintf(fid, format_spec, "",  "", "RMSE (mm)", ...
    statistics.EUC.RMSE_ve, statistics.AFF.RMSE_ve, statistics.TPS5.RMSE_ve, statistics.Kernel.RMSE_ve);
exportgraphics(hfig, [fdir, 'TOPACS_GPA_Visualization_4_356', '.pdf'], 'ContentType','vector');


fprintf(fid, "\\midrule \n");

trainingset_indicator = [5];    testset_indicator = [3, 4, 6];
[statistics, hfig] = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator);
fprintf(fid, format_spec, "",  "", "min (mm)", ...
    statistics.EUC.min, statistics.AFF.min, statistics.TPS5.min, statistics.Kernel.min);
fprintf(fid, format_spec, sprintf(class_str, trainingset_indicator),  sprintf(class_str, testset_indicator), "max (mm)", ...
    statistics.EUC.max, statistics.AFF.max, statistics.TPS5.max, statistics.Kernel.max);
fprintf(fid, format_spec, "",  "", "mean (mm)", ...
    statistics.EUC.mean, statistics.AFF.mean, statistics.TPS5.mean, statistics.Kernel.mean);
fprintf(fid, format_spec, "",  "", "RMSE (mm)", ...
    statistics.EUC.RMSE_ve, statistics.AFF.RMSE_ve, statistics.TPS5.RMSE_ve, statistics.Kernel.RMSE_ve);
exportgraphics(hfig, [fdir, 'TOPACS_GPA_Visualization_5_346', '.pdf'], 'ContentType','vector');


fprintf(fid, "\\midrule \n");

trainingset_indicator = [6];    testset_indicator = [3, 4, 5];
[statistics, hfig] = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator);
fprintf(fid, format_spec, "",  "", "min (mm)", ...
    statistics.EUC.min, statistics.AFF.min, statistics.TPS5.min, statistics.Kernel.min);
fprintf(fid, format_spec, sprintf(class_str, trainingset_indicator),  sprintf(class_str, testset_indicator), "max (mm)", ...
    statistics.EUC.max, statistics.AFF.max, statistics.TPS5.max, statistics.Kernel.max);
fprintf(fid, format_spec, "",  "", "mean (mm)", ...
    statistics.EUC.mean, statistics.AFF.mean, statistics.TPS5.mean, statistics.Kernel.mean);
fprintf(fid, format_spec, "",  "", "RMSE (mm)", ...
    statistics.EUC.RMSE_ve, statistics.AFF.RMSE_ve, statistics.TPS5.RMSE_ve, statistics.Kernel.RMSE_ve);
exportgraphics(hfig, [fdir, 'TOPACS_GPA_Visualization_6_345', '.pdf'], 'ContentType','vector');


fprintf(fid, "\\bottomrule \n");
fprintf(fid, "\\end{tabular}");
fclose (fid);


end



% the final case to save as an example
%trainingset_indicator = [3];    testset_indicator = [4, 5, 6];
%trainingset_indicator = [4];    testset_indicator = [4, 5, 6];
%trainingset_indicator = [6];    testset_indicator = [3, 4, 5, 6];
[statistics, hfig] = runTrainigTestTOPACS (DataCell, VisVecCell, [4], [5, 6], fdir);
exportgraphics(hfig, [fdir, 'TOPACS_GPA_Visualization_4_56', '.pdf'], 'ContentType','vector');



close all;

sampleChoice = 1;

tmp = readmatrix([datafiledir,'points', num2str(sampleChoice-1), '.csv']);
PointCloudsSample = tmp(:, 1:3)';
clear tmp
landmarks = DataCell{sampleChoice};
landmarks = landmarks(:, logical(VisVecCell{sampleChoice}));


hFig0 = figure("Name", 'TOPACS', 'Position', [100, 50, 400, 350]);
tfig0 =  tiledlayout(3, 2, 'TileSpacing','compact');

color1 = '#103232'
color2 = '#3742FA';
color3 = '#05C46B'
color4 = [0 1 1];
color5 = '#FC5C65'

% color2 = '#104FFF'
% color3 = '#2FD151'
% color4 = '#64C7B8'
% color5 = '#FF1038'
% 
% color5 = '#2FD151'
% color4 = '#64C7B8'
% color3 = '#FF1038'
% color2 = '#45CAFF'



ax2 = nexttile(tfig0, [3, 1]);
ax2.FontSize = 5;
ax2.LineWidth = 0.5;
hold on;
hpoints = scatter3(PointCloudsSample(1,:), PointCloudsSample(2,:), PointCloudsSample(3,:), 10, "Marker",".", "MarkerEdgeColor", color1, "MarkerFaceColor", color1);
corresp_points_cell = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, [3]);  corresp_points = corresp_points_cell{1};
c3 = scatter3(corresp_points(1,:), corresp_points(2,:), corresp_points(3,:), 4, "MarkerEdgeColor", color2, "MarkerFaceColor", color2);
corresp_points_cell = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, [4]);  corresp_points = corresp_points_cell{1};
c4 = scatter3(corresp_points(1,:), corresp_points(2,:), corresp_points(3,:), 4, "MarkerEdgeColor", color3, "MarkerFaceColor", color3);
corresp_points_cell = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, [5]);  corresp_points = corresp_points_cell{1};
c5 = scatter3(corresp_points(1,:), corresp_points(2,:), corresp_points(3,:), 4, "MarkerEdgeColor", color4, "MarkerFaceColor", color4);
corresp_points_cell = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, [6]);  corresp_points = corresp_points_cell{1};
c6 = scatter3(corresp_points(1,:), corresp_points(2,:), corresp_points(3,:), 4, "MarkerEdgeColor", color5, "MarkerFaceColor", color5);
hold off; view(3); axis equal tight; 
ax2.TickLength = [0.01, 0.01];
view(-170, 5)

lgd = legend([hpoints, c3, c4, c5, c6], {'point-cloud', '$\mathcal{C}_3$', '$\mathcal{C}_4$', '$\mathcal{C}_5$', '$\mathcal{C}_6$'}, 'FontSize', 10, 'Interpreter', 'latex', 'Orientation', 'vertical', 'box', 'off');
lgd.Layout.Tile = 2;
lgd.Layout.TileSpan = [2, 1];

pause(0.2);
exportgraphics(hFig0, [fdir, 'TOPACS_Sample_PointCloud_Visualization', '.pdf'], 'ContentType','vector');





%%% --- point clouds

trainingset_indicator = [6];    testset_indicator = [3, 4, 5];
[statistics, hfig, GPA_RIGID, GPA_AFFINE, GPA_TPS, GPA_KERNEL] = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator, fdir);


PointClouds = cell(1, 6);
for ii = 1 : 6
tmp = readmatrix([datafiledir,'points', num2str(ii-1), '.csv']);
PointClouds{ii} = tmp(:, 1:3)';
end




close all;


% GPA_RIGID
figure;
hold on;
for ii = 1 : 6
    transPts = GPA_RIGID.transformPoints(PointClouds{ii}, ii);
    hpoints = scatter3(transPts(1,:), transPts(2,:), transPts(3,:), 10, "Marker",".");
end
hold off;
view(3); axis equal tight


% GPA_AFFINE
figure;
hold on;
for ii = 1 : 6
    transPts = GPA_AFFINE.transformPoints(PointClouds{ii}, ii);
    hpoints = scatter3(transPts(1,:), transPts(2,:), transPts(3,:), 10, "Marker",".");
end
hold off;
view(3); axis equal tight


% GPA_TPS
figure;
hold on;
for ii = 1 : 6
    transPts = GPA_TPS.transformPoints(PointClouds{ii}, ii);
    hpoints = scatter3(transPts(1,:), transPts(2,:), transPts(3,:), 10, "Marker",".");
end
hold off;
view(3); axis equal tight


% GPA_KERNEL
figure;
hold on;
for ii = 1 : 6
    transPts = GPA_KERNEL.transformPoints(PointClouds{ii}, ii);
    hpoints = scatter3(transPts(1,:), transPts(2,:), transPts(3,:), 10, "Marker",".");
end
hold off;
view(3); axis equal tight





function [statistics, hfig, GPA_RIGID, GPA_AFFINE, GPA_TPS, GPA_KERNEL] = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator, fdir)

flagJointlyOptimizeLamdaRigidTransformation = true;

TPS_SmoothParam = 0.01;

Kernel_BandwithParam = 0.2;
Kernel_SmoothParam = 0.05;

figSize = [250, 300];


% We use correspondences that occur in 3, 4, 5 point-clouds to compute tranfsormations
[TrainingDataCell, TrainingVisVecCell] = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, trainingset_indicator);
[TestDataCell, TestVisVecCell] = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, testset_indicator);

% Rigid-GPA
GPA_RIGID = RigidGPA;
GPA_RIGID.addDataByArray(TrainingDataCell, TrainingVisVecCell);
GPA_RIGID.run();
GPA_RIGID


% Affine-GPA
GPA_AFFINE = DefGPA('AFFINE');
GPA_AFFINE.dimFeatureSpace = 0;
GPA_AFFINE.addDataByArray(TrainingDataCell, TrainingVisVecCell);
GPA_AFFINE.run (TPS_SmoothParam, flagJointlyOptimizeLamdaRigidTransformation);
GPA_AFFINE.rsdError();
GPA_AFFINE


% TPS-GPA. TPS(5)
GPA_TPS = DefGPA('TPS');
GPA_TPS.dimFeatureSpace =5^3;
GPA_TPS.addDataByArray(TrainingDataCell, TrainingVisVecCell);
GPA_TPS.run (TPS_SmoothParam, flagJointlyOptimizeLamdaRigidTransformation);
GPA_TPS.rsdError();
GPA_TPS


% Kernel-GPA
GPA_KERNEL = KernelGPA;
GPA_KERNEL.smoothParam = Kernel_SmoothParam; % 0.05

%fGPA_KERNEL.quantilePercentParam = Kernel_BandwithParam;
GPA_KERNEL.kernelScaleParam= Kernel_BandwithParam;  % 0.2

GPA_KERNEL.addDataByArray(TrainingDataCell, TrainingVisVecCell);
GPA_KERNEL.run (flagJointlyOptimizeLamdaRigidTransformation);
GPA_KERNEL

                

gaugeShape = GPA_RIGID.mShape;

tfontsize = 9;

close all;

hfig = figure("Position", [100, 100, 480,  600]);
tfig = tiledlayout(2, 2, "TileSpacing","compact");




nexttile;
% transform the point-clouds to the reference frame
rect_R = OrhoProcurestes (GPA_RIGID.mShape, gaugeShape);
for ii = 1 : length(TestDataCell)
    data = GPA_RIGID.transformPoints(TestDataCell{ii}, ii);
    transTestDataCell{ii} = rect_R * data;
end
set(gca, 'LineWidth', 0.5);
stdvar1 = VisualizePointCloudsVariation (transTestDataCell, TestVisVecCell);
view(3); view(86.5261,   2.9606); 
axis on; axis equal tight; box on; axis off;
pause(0.1);
%set(tfig1, 'Position', [100, 800, figSize(1), figSize(2)]);
set(gca, 'FontSize', 6);
title({'Rigid-GPA', ['max = ', num2str(stdvar1.max, '%.2f'), ' min = ', num2str(stdvar1.min, '%.2f'), ' mean = ', num2str(stdvar1.mean, '%.2f')], ['RMSE = ', num2str(stdvar1.RMSE_ve, '%.2f')]}, ...
    'Interpreter', 'Latex', 'FontSize', tfontsize);
% if exist("fdir", "var")
% exportgraphics(tfig1, [fdir, 'TOPACS_GPA_Visualization_Rigid', '.pdf'], 'ContentType','vector');
% end
pause(0.5);


xaxislim = xlim;
yaxislim = ylim;
zaxislim = zlim;

% papersize = tfig1.PaperSize;
% papersize(1) = 5.5;



nexttile;
% transform the point-clouds to the reference frame
rect_R = OrhoProcurestes (GPA_AFFINE.mShape, gaugeShape);
for ii = 1 : length(TestDataCell)
    data = GPA_AFFINE.transformPoints(TestDataCell{ii}, ii);
    transTestDataCell{ii} = rect_R * data;
end
stdvar2 = VisualizePointCloudsVariation (transTestDataCell, TestVisVecCell);
set(gca, 'LineWidth', 0.5);
view(3); view(86.5261,   2.9606);
axis on; axis equal tight; box on; axis off;
pause(0.1);
set(gca, 'XLim', xaxislim); set(gca, 'YLim', yaxislim); set(gca, 'ZLim', zaxislim);
%set(tfig2, 'Position', [400, 800, figSize(1), figSize(2)]);
set(gca, 'FontSize', 6);
title({'Affine-GPA', ['max = ', num2str(stdvar2.max, '%.2f'), ' min = ', num2str(stdvar2.min, '%.2f'), ' mean = ', num2str(stdvar2.mean, '%.2f')], ['RMSE = ', num2str(stdvar2.RMSE_ve, '%.2f')]}, ...
    'Interpreter', 'Latex', 'FontSize', tfontsize);
pause(0.5);




nexttile;
% transform the point-clouds to the reference frame
rect_R = OrhoProcurestes (GPA_TPS.mShape, gaugeShape);
for ii = 1 : length(TestDataCell)
    data = GPA_TPS.transformPoints(TestDataCell{ii}, ii);
    transTestDataCell{ii} = rect_R * data;
end
stdvar4 = VisualizePointCloudsVariation (transTestDataCell, TestVisVecCell);
set(gca, 'LineWidth', 0.5);
view(3); view(86.5261,   2.9606);
axis on; axis equal tight; box on; axis off;
pause(0.1);
set(gca, 'XLim', xaxislim); set(gca, 'YLim', yaxislim); set(gca, 'ZLim', zaxislim);
%set(tfig4, 'Position', [100, 100, figSize(1), figSize(2)]);
set(gca, 'FontSize', 6);
title({'TPS-GPA', ['max = ', num2str(stdvar4.max, '%.2f'), ' min = ', num2str(stdvar4.min, '%.2f'), ' mean = ', num2str(stdvar4.mean, '%.2f')], ['RMSE = ', num2str(stdvar4.RMSE_ve, '%.2f')]}, ...
    'Interpreter', 'Latex', 'FontSize', tfontsize);
pause(0.5);



nexttile;
% transform the point-clouds to the reference frame
rect_R = OrhoProcurestes (GPA_KERNEL.mShape, gaugeShape);
for ii = 1 : length(TestDataCell)
    data = GPA_KERNEL.transformPoints(TestDataCell{ii}, ii);
    transTestDataCell{ii} = rect_R * data;
end
set(gca, 'LineWidth', 0.5);
stdvar6 = VisualizePointCloudsVariation (transTestDataCell, TestVisVecCell);
view(3); view(86.5261,   2.9606);
axis on; axis equal tight; box on; axis off;
pause(0.1);
set(gca, 'XLim', xaxislim); set(gca, 'YLim', yaxislim); set(gca, 'ZLim', zaxislim);
set(gca, 'FontSize', 6);
title({'Kernel-GPA', ['max = ', num2str(stdvar6.max, '%.2f'), ' min = ', num2str(stdvar6.min, '%.2f'), ' mean = ', num2str(stdvar6.mean, '%.2f')], ['RMSE = ', num2str(stdvar6.RMSE_ve, '%.2f')]}, ...
    'Interpreter', 'Latex', 'FontSize', tfontsize);
pause(0.5);


cb = colorbar;
cb.Layout.Tile ="north";
pause(0.5);


statistics.EUC = stdvar1;
statistics.AFF = stdvar2;
statistics.TPS3 = 0;
statistics.TPS5 = stdvar4;
statistics.TPS7 = 0;
statistics.Kernel = stdvar6;





end






function signs = GetRectSigns (Shape1, ShapeGauge)

v = sum(Shape1, 2);
vg = sum(ShapeGauge, 2);

ptc = 2;

v =  Shape1(: ,ptc);
vg = ShapeGauge(: ,ptc);

signs = sign(v) ./ sign(vg);


end



% orthogonal Procrustes for  R * D1 = D2
function R = OrhoProcurestes (D1, D2)

meanD1 = mean(D1, 2);

meanD2 = mean(D2, 2);

M = (D1 - meanD1) * (D2 - meanD2)';

[U, ~, V] = svd(M);

R = V * U';

end

