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



datafiledir = '../../dataset/TOPACS/landmarks/';


class_str = "$\\mathcal{C}_{%d}$";


[DataCell, VisVecCell] = readDataset (datasetName);


fid = fopen ('../../PaperDraft/figure/exp_TOPACS_validation_table.tex', 'w');
format_spec = " %s  &  %s  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  \\textbf{%.2f} \\\\ \n";
fprintf(fid, "\\begin{tabular}{@{} cccccccc @{}} \n");
fprintf(fid, "\\toprule \n");
fprintf(fid, " training  &  test  &  EUC & AFF & TPS($3$) & TPS($5$) & TPS($7$)  &  Kernel  \\\\ \n");

fprintf(fid, "\\midrule \n");

trainingset_indicator = [3];    testset_indicator = [4, 5, 6];
statistics = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator);
fprintf(fid, format_spec, sprintf(class_str, trainingset_indicator),  sprintf(class_str, testset_indicator), ...
    statistics.EUC, statistics.AFF, statistics.TPS3, statistics.TPS5, statistics.TPS7, statistics.Kernel);

fprintf(fid, "\\midrule \n");

trainingset_indicator = [4];    testset_indicator = [3, 5, 6];
statistics = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator);
fprintf(fid, format_spec, sprintf(class_str, trainingset_indicator),  sprintf(class_str, testset_indicator), ...
    statistics.EUC, statistics.AFF, statistics.TPS3, statistics.TPS5, statistics.TPS7, statistics.Kernel);

fprintf(fid, "\\midrule \n");

trainingset_indicator = [5];    testset_indicator = [3, 4, 6];
statistics = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator);
fprintf(fid, format_spec, sprintf(class_str, trainingset_indicator),  sprintf(class_str, testset_indicator), ...
    statistics.EUC, statistics.AFF, statistics.TPS3, statistics.TPS5, statistics.TPS7, statistics.Kernel);

fprintf(fid, "\\midrule \n");

trainingset_indicator = [6];    testset_indicator = [3, 4, 5];
statistics = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator);
fprintf(fid, format_spec, sprintf(class_str, trainingset_indicator),  sprintf(class_str, testset_indicator), ...
    statistics.EUC, statistics.AFF, statistics.TPS3, statistics.TPS5, statistics.TPS7, statistics.Kernel);

fprintf(fid, "\\bottomrule \n");
fprintf(fid, "\\end{tabular}");
fclose (fid);


% the final case to save as an example
trainingset_indicator = [4];    testset_indicator = [5, 6];
statistics = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator, fdir);





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




function statistics = runTrainigTestTOPACS (DataCell, VisVecCell, trainingset_indicator, testset_indicator, fdir)

flagJointlyOptimizeLamdaRigidTransformation = true;

TPS_SmoothParam = 0.01;

Kernel_Bandwith_Quantile = 0.2;
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
GPA_TPS3 = DefGPA('TPS');
GPA_TPS3.dimFeatureSpace =3^3;
GPA_TPS3.addDataByArray(TrainingDataCell, TrainingVisVecCell);
GPA_TPS3.run (TPS_SmoothParam, flagJointlyOptimizeLamdaRigidTransformation);
GPA_TPS3.rsdError();
GPA_TPS3


% TPS-GPA. TPS(5)
GPA_TPS5 = DefGPA('TPS');
GPA_TPS5.dimFeatureSpace =5^3;
GPA_TPS5.addDataByArray(TrainingDataCell, TrainingVisVecCell);
GPA_TPS5.run (TPS_SmoothParam, flagJointlyOptimizeLamdaRigidTransformation);
GPA_TPS5.rsdError();
GPA_TPS5


% TPS-GPA. TPS(7)
GPA_TPS7 = DefGPA('TPS');
GPA_TPS7.dimFeatureSpace =7^3;
GPA_TPS7.addDataByArray(TrainingDataCell, TrainingVisVecCell);
GPA_TPS7.run (TPS_SmoothParam, flagJointlyOptimizeLamdaRigidTransformation);
GPA_TPS7.rsdError();
GPA_TPS7


% Kernel-GPA
GPA_KERNEL = KernelGPA;
GPA_KERNEL.smoothParam = Kernel_SmoothParam;
GPA_KERNEL.quantilePercentParam= Kernel_Bandwith_Quantile;
GPA_KERNEL.addDataByArray(TrainingDataCell, TrainingVisVecCell);
GPA_KERNEL.run (flagJointlyOptimizeLamdaRigidTransformation);
GPA_KERNEL

                

gaugeShape = GPA_RIGID.mShape;

tfontsize = 7;

close all;

hfig = figure("Position", [100, 100, 480,  500]);
tfig = tiledlayout(2, 3, "TileSpacing","compact");

nexttile;
%tfig1 = figure(1);
set(gca, 'LineWidth', 0.5);
% transform the point-clouds to the reference frame
rect_R = OrhoProcurestes (GPA_RIGID.mShape, gaugeShape);
for ii = 1 : length(TestDataCell)
    data = GPA_RIGID.transformPoints(TestDataCell{ii}, ii);
    transTestDataCell{ii} = rect_R * data;
end
stdvar1 = VisualizePointCloudsVariation (transTestDataCell, TestVisVecCell);
view(-91.4467,   -1.0943); axis on; axis equal tight; box on; axis off;
pause(0.1);
%set(tfig1, 'Position', [100, 800, figSize(1), figSize(2)]);
set(gca, 'FontSize', 6);
title({'EUC', ['std-dev= ', num2str(stdvar1, '%.2f')]}, 'Interpreter', 'Latex', 'FontSize', tfontsize);
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
%tfig2 = figure(2);
set(gca, 'LineWidth', 0.5);
% transform the point-clouds to the reference frame
rect_R = OrhoProcurestes (GPA_AFFINE.mShape, gaugeShape);
for ii = 1 : length(TestDataCell)
    data = GPA_AFFINE.transformPoints(TestDataCell{ii}, ii);
    transTestDataCell{ii} = rect_R * data;
end
stdvar2 = VisualizePointCloudsVariation (transTestDataCell, TestVisVecCell);
view(-91.4467,   -1.0943); axis on; axis equal tight; box on; axis off;
pause(0.1);
set(gca, 'XLim', xaxislim); set(gca, 'YLim', yaxislim); set(gca, 'ZLim', zaxislim);
%set(tfig2, 'Position', [400, 800, figSize(1), figSize(2)]);
set(gca, 'FontSize', 6);
title({'AFF', ['std-dev= ', num2str(stdvar2, '%.2f')]}, 'Interpreter', 'Latex', 'FontSize', tfontsize);
% if exist("fdir", "var")
% exportgraphics(tfig2, [fdir, 'TOPACS_GPA_Visualization_Affine', '.pdf'], 'ContentType','vector');
% end
pause(0.5);


nexttile;
%tfig3 = figure(3);
set(gca, 'LineWidth', 0.5);
% transform the point-clouds to the reference frame
rect_R = OrhoProcurestes (GPA_TPS3.mShape, gaugeShape);
for ii = 1 : length(TestDataCell)
    data = GPA_TPS3.transformPoints(TestDataCell{ii}, ii);
    transTestDataCell{ii} = rect_R * data;
end
stdvar3 = VisualizePointCloudsVariation (transTestDataCell, TestVisVecCell);
view(-91.4467,   -1.0943); axis on; axis equal tight; box on; axis off;
pause(0.1);
set(gca, 'XLim', xaxislim); set(gca, 'YLim', yaxislim); set(gca, 'ZLim', zaxislim);
%set(tfig3, 'Position', [700, 800, figSize(1), figSize(2)]);
set(gca, 'FontSize', 6);
title({'TPS(3)', ['std-dev= ', num2str(stdvar3, '%.2f')]},'Interpreter', 'Latex', 'FontSize', tfontsize);
% if exist("fdir", "var")
% exportgraphics(tfig3, [fdir, 'TOPACS_GPA_Visualization_TPS3', '.pdf'], 'ContentType','vector');
% end
pause(0.5);


nexttile;
%tfig4 = figure(4);
set(gca, 'LineWidth', 0.5);
% transform the point-clouds to the reference frame
rect_R = OrhoProcurestes (GPA_TPS5.mShape, gaugeShape);
for ii = 1 : length(TestDataCell)
    data = GPA_TPS5.transformPoints(TestDataCell{ii}, ii);
    transTestDataCell{ii} = rect_R * data;
end
stdvar4 = VisualizePointCloudsVariation (transTestDataCell, TestVisVecCell);
view(-91.4467,   -1.0943); axis on; axis equal tight; box on; axis off;
pause(0.1);
set(gca, 'XLim', xaxislim); set(gca, 'YLim', yaxislim); set(gca, 'ZLim', zaxislim);
%set(tfig4, 'Position', [100, 100, figSize(1), figSize(2)]);
set(gca, 'FontSize', 6);
title({'TPS(5)', ['std-dev= ', num2str(stdvar4, '%.2f')]}, 'Interpreter', 'Latex', 'FontSize', tfontsize);
% if exist("fdir", "var")
% exportgraphics(tfig4, [fdir, 'TOPACS_GPA_Visualization_TPS5', '.pdf'], 'ContentType','vector');
% end
pause(0.5);


nexttile;
%tfig5 = figure(5);
set(gca, 'LineWidth', 0.5);
% transform the point-clouds to the reference frame
rect_R = OrhoProcurestes (GPA_TPS7.mShape, gaugeShape);
for ii = 1 : length(TestDataCell)
    data = GPA_TPS7.transformPoints(TestDataCell{ii}, ii);
    transTestDataCell{ii} = rect_R * data;
end
stdvar5 = VisualizePointCloudsVariation (transTestDataCell, TestVisVecCell);
view(-91.4467,   -1.0943); axis on; axis equal tight; box on; axis off;
pause(0.1);
set(gca, 'XLim', xaxislim); set(gca, 'YLim', yaxislim); set(gca, 'ZLim', zaxislim);
%set(tfig5, 'Position', [400, 100, figSize(1), figSize(2)]);
set(gca, 'FontSize', 6);
title({'TPS(7)', ['std-dev= ', num2str(stdvar5, '%.2f')]}, 'Interpreter', 'Latex', 'FontSize', tfontsize);
% if exist("fdir", "var")
% exportgraphics(tfig5, [fdir, 'TOPACS_GPA_Visualization_TPS7', '.pdf'], 'ContentType','vector');
% end
pause(0.5);


nexttile;
%tfig6 = figure(6);
set(gca, 'LineWidth', 0.5);
% transform the point-clouds to the reference frame
rect_R = OrhoProcurestes (GPA_KERNEL.mShape, gaugeShape);
for ii = 1 : length(TestDataCell)
    data = GPA_KERNEL.transformPoints(TestDataCell{ii}, ii);
    transTestDataCell{ii} = rect_R * data;
end
stdvar6 = VisualizePointCloudsVariation (transTestDataCell, TestVisVecCell);
view(-91.4467,   -1.0943); axis on; axis equal tight; box on; axis off;
pause(0.1);
set(gca, 'XLim', xaxislim); set(gca, 'YLim', yaxislim); set(gca, 'ZLim', zaxislim);
%set(tfig6, 'Position', [700, 100, figSize(1), figSize(2)]);
set(gca, 'FontSize', 6);
title({'Kernel', ['std-dev= ', num2str(stdvar6, '%.2f')]}, 'Interpreter', 'Latex', 'FontSize', tfontsize);
% if exist("fdir", "var")
% exportgraphics(tfig6, [fdir, 'TOPACS_GPA_Visualization_Kernel', '.pdf'], 'ContentType','vector');
% end
pause(0.5);



if exist("fdir", "var")
exportgraphics(hfig, [fdir, 'TOPACS_GPA_Visualization', '.pdf'], 'ContentType','vector');
end




statistics.EUC = stdvar1;
statistics.AFF = stdvar2;
statistics.TPS3 = stdvar3;  
statistics.TPS5 = stdvar4;
statistics.TPS7 = stdvar5;
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

