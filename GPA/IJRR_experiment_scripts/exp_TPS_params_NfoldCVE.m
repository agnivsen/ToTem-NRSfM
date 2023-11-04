clear all
close all
clc


fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end

ffdir = '../../dataset';

% include the files in DefGPA directory
addpath('../DefGPA', '../');


datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug', 'TOPACS'}


dataOptions = [1  2  3  4, 6, 7];

Nfold = 20;    
flagJointlyOptimizeLamdaRigidTransformation = true;


dataOptions = 1

for dtc = dataOptions

datasetName = datasetCell{dtc};


[DataCell, VisVecCell] = readDataset (datasetName, ffdir);


if strcmp(datasetName, 'TOPACS')
    % We use correspondences that occur in 3, 4, 5 point-clouds to compute tranfsormations
    [DataCell, VisVecCell] = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, [4]);
end

if size(DataCell{1}, 1) == 3
    TPSParamRange = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1];
end
if size(DataCell{1}, 1) == 2
    TPSParamRange = [0.01, 0.1, 1:1:9, 10:10:300];
end
 

fprintf('\nDefPGA - TPS(5)\n');
GPA_TPS5 = DefGPA('TPS');
GPA_TPS5.dimFeatureSpace = 5^3;
GPA_TPS5.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
GPA_TPS5.addDataByArray(DataCell, VisVecCell);
GPA_TPS5.verbose = true;
arrayTrainingStatisticsNfold5 = GPA_TPS5.OptimizeSmoothParam (TPSParamRange, Nfold);



TPS5_PARAM_Nfold_CELL{dtc} = arrayTrainingStatisticsNfold5;


save(['TPS_PARAM_',datasetName, '_' num2str(Nfold), 'Fold.mat'], 'TPS5_PARAM_Nfold_CELL');



close all;


arrayTrainingStatisticsNfold5 = TPS5_PARAM_Nfold_CELL{dtc};


% visualization
figsize = [450, 240];
linewidth = 1.0;

hFig = figure(1000);


color1 = '#EE0000'
color2 = '#00EE00'
color3 = '#0000EE'



% color1 ='c';
% color2 = '#50F';


optSmoothParam = 0.01;
params = arrayTrainingStatisticsNfold5.SmoothParam;

ax1 = gca;
ax1.FontSize = 7;
hold on;
h1 = plot(params, arrayTrainingStatisticsNfold5.RMSE_ref, 'LineWidth', linewidth, 'Color', color1, 'LineStyle', '-', 'Marker', '.');
h2 = plot(params, arrayTrainingStatisticsNfold5.RMSE_dat, 'LineWidth', linewidth, 'Color', color2, 'LineStyle', '-.', 'Marker', 'square');
h3 = plot(params, arrayTrainingStatisticsNfold5.CVE, 'LineWidth', linewidth, 'Color', color3, 'LineStyle', '-.', 'Marker', 's');
hold off;
box on;
ahx = xline (optSmoothParam, '--r', {['\mu = ', num2str(optSmoothParam)]}, 'LineWidth', linewidth, 'FontSize', 8, 'LabelVerticalAlignment', 'middle', "LabelHorizontalAlignment", 'center');
set(ax1, 'xscale','log');
legend([h1, h2, h3], {'TPS(5)-RMSE\_r', 'TPS(5)-RMSE\_d', sprintf('TPS(5)-CVE-%dfold', Nfold)}, 'FontSize', 7, 'box', 'off',  'Interpreter', "latex");
ylabel('RMSE/CVE',  'FontSize', 8, 'Interpreter', "latex");
xlabel('regularization strength $\mu$', 'FontSize', 10, 'Interpreter', 'latex');


hFig.Position = [0, 800, figsize(1), figsize(2)];
exportgraphics(hFig, [fdir, 'GPA_TPS_Param_Visualization_', datasetName '.pdf'], 'ContentType','vector');


end





