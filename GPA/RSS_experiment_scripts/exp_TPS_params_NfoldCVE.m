clear all
close all
clc


fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end

% include the files in DefGPA directory
addpath('../DefGPA', '../');


datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug', 'TOPACS'}


dataOptions = [1  2  3  4, 6, 7];

Nfold = 20;    
flagJointlyOptimizeLamdaRigidTransformation = true;


dataOptions = 7

for dtc = dataOptions

datasetName = datasetCell{dtc};


[DataCell, VisVecCell] = readDataset (datasetName);


if strcmp(datasetName, 'TOPACS')
    % We use correspondences that occur in 3, 4, 5 point-clouds to compute tranfsormations
    [DataCell, VisVecCell] = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, [4]);
end

if size(DataCell{1}, 1) == 3
    TPSParamRange = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1:0.1:1];
end
if size(DataCell{1}, 1) == 2
    TPSParamRange = [0.01, 0.1, 1:1:9, 10:10:300];
end
 
% GPA_TPS3 = DefGPA('TPS');
% GPA_TPS3.dimFeatureSpace = 3^3;
% GPA_TPS3.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
% GPA_TPS3.addDataByArray(DataCell, VisVecCell);
% GPA_TPS3.verbose = false;
% arrayTrainingStatistics3 = GPA_TPS3.OptimizeSmoothParam (TPSParamRange, Nfold);

fprintf('\nDefPGA - TPS(5)\n');
GPA_TPS5 = DefGPA('TPS');
GPA_TPS5.dimFeatureSpace = 5^3;
GPA_TPS5.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
GPA_TPS5.addDataByArray(DataCell, VisVecCell);
GPA_TPS5.verbose = false;
arrayTrainingStatisticsNfold5 = GPA_TPS5.OptimizeSmoothParam (TPSParamRange, Nfold);
arrayTrainingStatisticsLeaveOne5 = GPA_TPS5.OptimizeSmoothParam (TPSParamRange, inf);


fprintf('\nDefPGA - TPS(7)\n');
GPA_TPS7 = DefGPA('TPS');
GPA_TPS7.dimFeatureSpace = 7^3;
GPA_TPS7.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
GPA_TPS7.addDataByArray(DataCell, VisVecCell);
GPA_TPS7.verbose = false;
arrayTrainingStatisticsNfold7 = GPA_TPS7.OptimizeSmoothParam (TPSParamRange, Nfold);
arrayTrainingStatisticsLeaveOne7 = GPA_TPS7.OptimizeSmoothParam (TPSParamRange, inf);



TPS5_PARAM_Nfold_CELL{dtc} = arrayTrainingStatisticsNfold5;
TPS7_PARAM_Nfold_CELL{dtc} = arrayTrainingStatisticsNfold7;

TPS5_PARAM_LeaveOne_CELL{dtc} = arrayTrainingStatisticsLeaveOne5;
TPS7_PARAM_LeaveOne_CELL{dtc} = arrayTrainingStatisticsLeaveOne7;


save(['TPS_PARAM_', num2str(Nfold), 'Fold.mat'], 'TPS5_PARAM_Nfold_CELL', 'TPS7_PARAM_Nfold_CELL', 'TPS5_PARAM_LeaveOne_CELL', 'TPS7_PARAM_LeaveOne_CELL');



close all;

dtc = 7;

arrayTrainingStatisticsNfold5 = TPS5_PARAM_Nfold_CELL{dtc};
arrayTrainingStatisticsNfold7 = TPS7_PARAM_Nfold_CELL{dtc};

arrayTrainingStatisticsLeaveOne5 = TPS5_PARAM_LeaveOne_CELL{dtc};
arrayTrainingStatisticsLeaveOne7 = TPS7_PARAM_LeaveOne_CELL{dtc};

% visualization
figsize = [450, 240];
linewidth = 1.0;

hFig = figure(1000);
htfig = tiledlayout(1, 2);
htfig.TileSpacing = 'compact';

c1 = '#4071EC'
c2 = '#18C1FF'


color1 = c1
color11 = c1
color2 = c2
color22 = c2


% color1 ='c';
% color2 = '#50F';


optSmoothParam = 0.01;
params = arrayTrainingStatisticsNfold5.SmoothParam;

ax1 = nexttile(htfig);
ax1.FontSize = 7;
hold on;
h1 = plot(params, arrayTrainingStatisticsNfold5.RMSE_ref, 'LineWidth', linewidth, 'Color', color1, 'LineStyle', '-', 'Marker', '.');
h2 = plot(params, arrayTrainingStatisticsNfold5.RMSE_dat, 'LineWidth', linewidth, 'Color', color11, 'LineStyle', '-.', 'Marker', 'square');
h3 = plot(params, arrayTrainingStatisticsNfold7.RMSE_ref, 'LineWidth', linewidth, 'Color', color2, 'LineStyle', '-', 'Marker', '.');
h4 = plot(params, arrayTrainingStatisticsNfold7.RMSE_dat, 'LineWidth', linewidth, 'Color', color22, 'LineStyle', '-.', 'Marker', 'square');
hold off;
box on;
ahx = xline (optSmoothParam, '--r', {['\mu = ', num2str(optSmoothParam)]}, 'LineWidth', linewidth, 'FontSize', 8, 'LabelVerticalAlignment', 'middle', "LabelHorizontalAlignment", 'center');
set(ax1, 'xscale','log');
legend([h1, h2, h3, h4], {'TPS(5)-RMSE\_r', 'TPS(5)-RMSE\_d', 'TPS(7)-RMSE\_r', 'TPS(7)-RMSE\_d'}, 'FontSize', 7, 'box', 'off',  'Interpreter', "latex");
ylabel('RMSE',  'FontSize', 8, 'Interpreter', "latex");


ax2 = nexttile(htfig);
ax2.FontSize = 7;
hold on;
h1 = plot(params, arrayTrainingStatisticsLeaveOne5.CVE, 'LineWidth', linewidth, 'Color', color1, 'LineStyle', '-', 'Marker', '.');
h2 = plot(params, arrayTrainingStatisticsNfold5.CVE, 'LineWidth', linewidth, 'Color', color1, 'LineStyle', '-.', 'Marker', 's');
h3 = plot(params, arrayTrainingStatisticsLeaveOne7.CVE, 'LineWidth', linewidth, 'Color', color2, 'LineStyle', '-', 'Marker', '.');
h4 = plot(params, arrayTrainingStatisticsNfold7.CVE, 'LineWidth', linewidth, 'Color', color2, 'LineStyle', '-.', 'Marker', 's');
hold off;
box on;
ahx = xline (optSmoothParam, '--r', {['\mu = ', num2str(optSmoothParam)]}, 'LineWidth', linewidth, 'FontSize', 8, 'LabelVerticalAlignment', 'middle', "LabelHorizontalAlignment", 'center');
set(ax2, 'xscale','log');
legend([h1, h2, h3, h4], {'TPS(5)-CVE-leave1',sprintf('TPS(5)-CVE-%dfold', Nfold), 'TPS(7)-CVE-leave1', sprintf('TPS(7)-CVE-%dfold', Nfold)}, 'FontSize', 7, 'box', 'off',  'Interpreter', "latex");
ylabel('CVE',  'Interpreter', 'latex', 'FontSize', 8);

xlabel(htfig, 'regularization strength $\mu$', 'FontSize', 10, 'Interpreter', 'latex');

hFig.Position = [0, 800, figsize(1), figsize(2)];
exportgraphics(hFig, [fdir, 'GPA_TPS_Param_Visualization_', datasetName '.pdf'], 'ContentType','vector');


end





