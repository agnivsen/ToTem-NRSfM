clear all
close all
clc


addpath('../DefGPA', '../');
addpath('../SGPAv1.0', '../SGPAv1.0/util');


figformat = 'pdf';


%datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug', 'TOPACS'}

datasetCell = {'LiverKanzhi', 'FacialExpressionKanzhi', 'TOPACS'}


fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end

dataOptions = [1  2  3];

Nfold = 20;


%smoothParams = [0.0001:0.0001:0.001, 0.002:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1, 2:1:10, 20:10:100];
%smoothParams = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100];


searchSmoothParams = true;

if (searchSmoothParams)
    smoothParams = [0.0001:0.0001:0.001, 0.002:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1];
    %smoothParams = [0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5]
    kernelScaleParams = [0.05 : 0.05 : 1];
else
    smoothParams = [0.001, 0.01, 0.1, 1];
    kernelScaleParams = [0.05 : 0.05 : 1];
end


for dtc = dataOptions

    ffdir = '../../dataset';
    datasetName = datasetCell{dtc};

    if strcmp(datasetName, 'LiverKanzhi')  
        load([ffdir, '/LiverKanzhi/liver_raw_dataset.mat']);
    end

    if strcmp(datasetName, 'FacialExpressionKanzhi')        
        load([ffdir, '/FacialExpressionKanzhi/face_raw_dataset.mat']);
    end

    if strcmp(datasetName, 'TOPACS')        
        [DataCell, VisVecCell] = readDataset (datasetName, ffdir);        
        % We use correspondences that occur in 3, 4, 5 point-clouds to compute tranfsormations
        [DataCell, VisVecCell] = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, [3]);
    end


    GPA_KERNEL = KernelGPA;
    GPA_KERNEL.verbose = false;
    GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
    GPA_KERNEL.flagJointlyOptimizeLamdaRigidTransformation = true;

    arrayTrainingStatistics = GPA_KERNEL.OptimizeHyperParams(smoothParams,  kernelScaleParams, Nfold);
    KERNEL_PARAM_STATISTICS_CELL {dtc} = arrayTrainingStatistics;

end




%%% ----- plot -----
%%% Data from Kanzhi ---
%%% Run cross-valiation for
% --- Liver
% --- Face
% --- TOPACS
% by setting:
%smoothParams = [0.0001:0.0001:0.001, 0.002:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1, 2:1:10, 20:10:100];
%kernelScaleParams = [0.25, 0.5, 0.75];
%
% Make sure you can generate three figures (one for each data) using the following script
% Then send me the structure KERNEL_PARAM_STATISTICS_CELL


close all;


%load('KERNEL_PARAM_10Fold200.mat');

for dtc = dataOptions

arrayTrainingStatistics = KERNEL_PARAM_STATISTICS_CELL{dtc};

v_mu = arrayTrainingStatistics.SmoothParam(:, 1); % column vector
v_p = arrayTrainingStatistics.KernelScaleParam(1, :); % row vector
v_CVE = arrayTrainingStatistics.CVE;

linestyles = {"-", "--", "-.", ":"};

if(searchSmoothParams)

    handles = [];
    lgdNames = {};
    hfig = figure('Position', [600, 600, 300, 200]);
    hold on;
    for kk = 1 : length(v_p)
        p = v_p(kk);
        CVE = v_CVE(:, kk);
        h = plot(v_mu, CVE, 'LineWidth', 1.5, "LineStyle", linestyles{kk});
        handles = [handles, h];
        lgdNames = [lgdNames, sprintf("$p$ = %.2f", p) ];
    end
    hold off;
    set(gca, "XScale", "log");
    lgd = legend(handles, lgdNames, 'FontSize', 14,'Orientation', 'Vertical', 'box', 'off', ...
            "Location", "northwest", "FontName", "Arial", "Interpreter", "Latex");
    xlabel("$\mu$", 'Fontsize', 14, "Interpreter", "Latex");    
    ylabel("\textrm{CVE}", 'Fontsize', 14, "Interpreter", "Latex");

    exportgraphics(hfig, [fdir, datasetName, '_mu.pdf'], 'ContentType','vector');
else

    handles = [];
    lgdNames = {};
    hfig = figure('Position', [600, 600, 300, 200]);
    hold on;
    for kk = 1 : length(v_mu)
        mu = v_mu(kk);
        CVE = v_CVE(kk, :);
        h = plot(v_p, CVE, 'LineWidth', 1.5, "LineStyle", linestyles{kk});
        handles = [handles, h];
        lgdNames = [lgdNames, sprintf("$\\mu$ = %.3f", mu) ];
    end
    hold off;
    set(gca, "XScale", "linear");
    lgd = legend(handles, lgdNames, 'FontSize', 10,'Orientation', 'Vertical', 'box', 'on', ...
            "Location", "northwest", "FontName", "Arial", "Interpreter", "Latex");
    xlabel("$p$", 'Fontsize', 14, "Interpreter", "Latex");    
    ylabel("\textrm{CVE}", 'Fontsize', 14, "Interpreter", "Latex");

    exportgraphics(hfig, [fdir, datasetName, '_p.pdf'], 'ContentType','vector');

end

end




% function cve = runKernelGPA_Nfold (DataCell, VisVecCell, Nfold)
% 
% GPA_KERNEL = KernelGPA;
% GPA_KERNEL.verbose = false;
% GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
% GPA_KERNEL.flagJointlyOptimizeLamdaRigidTransformation = true;
% 
% gaugeShape = GPA_KERNEL.mShape;
% 
% numPointClouds = GPA_KERNEL.numPointClouds;
% numPoints = GPA_KERNEL.numPoints;
% 
% 
% predictedPointClouds = cell(1, numPointClouds);
% for ii = 1 : this.numPointClouds
%     predictedPointClouds{ii} = 0 * gaugeShape;
% end
% N = min(numPoints, N);
% for ii = 1 : N
%     vtmp = (ii : N : numPoints);
%     testingPointsFlag = false(1, numPoints);
%     testingPointsFlag(vtmp) = true;
%     trainingPointsFlag = ~logical(testingPointsFlag);
%     % compute the transformed test points
%     [transform_points_cell, transform_vis_cell, refShape] = runTrainingTestData(this, trainingPointsFlag, testingPointsFlag);
%     % correct gauge and fill the predicted points
%     [sR, t] = this.EuclideanProcrustes (refShape, gaugeShape(:, trainingPointsFlag));
%     for kk = 1 : this.numPointClouds
%         predictedPointClouds{kk}(:, testingPointsFlag) = (sR * transform_points_cell{kk} + t) .* transform_vis_cell{kk};
%     end
% end
% 
% 
% end


