clear all
close all
clc


addpath('../DefGPA', '../');
addpath('../SGPAv1.0', '../SGPAv1.0/util');


figformat = 'pdf';


datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug', 'TOPACS'}

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end

dataOptions = [1  2  3  4, 6, 7];

Nfold = 20;
flagJointlyOptimizeLamdaRigidTransformation = true;


smoothParams = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100];
quantilePercentParams = [0.25, 0.5, 0.75];


for dtc = dataOptions

    ffdir = '../../dataset';

    datasetName = datasetCell{dtc}
    [DataCell, VisVecCell] = readDataset (datasetName, ffdir);
    
    if strcmp(datasetName, 'TOPACS')
        % We use correspondences that occur in 3, 4, 5 point-clouds to compute tranfsormations
        [DataCell, VisVecCell] = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, [4]);
    end

    GPA_KERNEL = KernelGPA;
    GPA_KERNEL.verbose = false;
    GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
    GPA_KERNEL.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
    
    arrayTrainingStatistics = GPA_KERNEL.OptimizeHyperParams(smoothParams,  quantilePercentParams, Nfold);
    KERNEL_PARAM_STATISTICS_CELL {dtc} = arrayTrainingStatistics;

end


close all
for dtc = dataOptions

    datasetName = datasetCell{dtc};
    arrayTrainingStatistics = KERNEL_PARAM_STATISTICS_CELL {dtc};

    qZt = arrayTrainingStatistics.CVE';
    qX = arrayTrainingStatistics.SmoothParam';
    qY = arrayTrainingStatistics.QuantilePercentParam';
    qRMSEt = arrayTrainingStatistics.RMSE_ref';

    hfig = figure('Name', datasetName, 'Position', [600, 600, 300, 300]); clf; hold on;    
    imagesc(qZt);
    set(gca, 'XTick',1:length(smoothParams), 'XTickLabel', smoothParams, 'XTickLabelRotation',30, 'FontSize', 8);
    set(gca, 'YTick',1:length(quantilePercentParams), 'YTickLabel', quantilePercentParams, 'YTickLabelRotation', 0, 'FontSize', 8);
    xlabel('$\mu$', 'Interpreter', 'latex');
    ylabel('$p$-quantile', 'Interpreter', 'latex');
    colormap(hot)
    colorbar
    hold off;
    title(datasetName);
    axis equal tight;
    
    exportgraphics(hfig, [fdir, 'KernelGPA_Params_', datasetName '.pdf'],'ContentType','vector');
    pause(0.1);

end

save(['KERNEL_PARAM_', num2str(Nfold), 'Fold.mat'], 'KERNEL_PARAM_STATISTICS_CELL');



close all
hfig = figure('Position', [200, 200, 600, 700]);
tfig = tiledlayout(ceil(length(dataOptions)/2), 2,"TileSpacing","tight");
for dtc = dataOptions
    
    datasetName = datasetCell{dtc};
    arrayTrainingStatistics = KERNEL_PARAM_STATISTICS_CELL {dtc};

    qZ = arrayTrainingStatistics.CVE;
    qX = arrayTrainingStatistics.SmoothParam;
    qY = arrayTrainingStatistics.QuantilePercentParam;
    qRMSE = arrayTrainingStatistics.RMSE_ref;

    ax = nexttile;
    ax.TickDir = 'out';
    ax.FontSize = 6;

    [X, Y] = meshgrid(1:length(smoothParams), 1:length(quantilePercentParams));

    mesh(X', Y', 0*qZ, qZ); hold on;
    colormap("turbo")
    scatter3(X(:), Y(:), 0*qZ(:), 50, qZ(:), 'filled'); hold off;
    colormap("turbo")
    view(0, 90);

    set(gca, 'XTick',1:length(smoothParams), 'XTickLabel', smoothParams, 'XTickLabelRotation',90, 'FontSize', 8);
    set(gca, 'YTick',1:length(quantilePercentParams), 'YTickLabel', quantilePercentParams, 'YTickLabelRotation', 0, 'FontSize', 8);


    colormap("turbo")
    cb = colorbar;
    cb.FontSize = 6;

    hold off;
    title(datasetName, 'FontWeight','normal', 'FontSize', 10);
    axis tight;
    grid off;
     
end
xlabel(tfig, '\textrm{kernel smooth parameter} $\mu$', 'Fontsize', 14, "Interpreter", "Latex");    
ylabel(tfig, '\textrm{kernel bandwith parameter} $p$-\textrm{quantile}', 'Fontsize', 14, "Interpreter", "Latex");
pause(0.1);
exportgraphics(hfig, [fdir, 'KernelGPA_Params', '.pdf'],'ContentType','vector');

