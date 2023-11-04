clear all
close all
clc

ffdir = '../../dataset';

addpath('../DefGPA', '../');
addpath('../SGPAv1.0', '../SGPAv1.0/util');



% the parameters to test

TPSsmoothParams = [1e-6, 1e-5, 1e-4, 0.001:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1];


smoothParams = [0.0001:0.0001:0.001, 0.001:0.001:0.01, 0.02:0.01:0.1, 0.2:0.1:1, 2:1:10, 20:10:100];
%kernelParams = [0.001, 0.01, 0.05, 0.1, 0.2];
kernelParams = [0.05 : 0.05 : 0.5];

datasetCell = {'LiverKanzhi', 'FacialExpressionKanzhi', 'TOPACS'}


kernelParamsStr = {};
for ii = 1 : length(kernelParams)
    kernelParamsStr{ii} = num2str(kernelParams(ii));
end



for dtc = [2, 3, 1]

    ffdir = '../../dataset';
    datasetName = datasetCell{dtc};

    if strcmp(datasetName, 'LiverKanzhi')
        load([ffdir, '/LiverKanzhi/liver_raw_dataset.mat']);
        TestPointsDataCell = LostDataCell;
        for k = 1 : length(LostDataCell)
            TestPointsVisVecCell{k} = ones(1, size(LostDataCell{k}, 2));
        end
    end

    if strcmp(datasetName, 'FacialExpressionKanzhi')
        load([ffdir, '/FacialExpressionKanzhi/face_raw_dataset.mat']);
        TestPointsDataCell = LostDataCell;
        for k = 1 : length(LostDataCell)
            TestPointsVisVecCell{k} = ones(1, size(LostDataCell{k}, 2));
        end        
    end


    if strcmp(datasetName, 'TOPACS')
        [DataCell, VisVecCell] = readDataset(datasetName, ffdir);
        [TestPointsDataCell, TestPointsVisVecCell] = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, [4, 5, 6]);
        [DataCell, VisVecCell] = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, [3]);
    end


    % test smooth parameters for Kernel
    Kernel_Point_Std_Dev_Matrix = [];

    for ii = 1 : length(kernelParams)
        for jj = 1 : length(smoothParams)
            [points_std_dev, RMSE_VE] = runTrainingTest_KernelGPA (DataCell, VisVecCell, TestPointsDataCell, TestPointsVisVecCell, ...
                smoothParams(jj), kernelParams(ii));
            Kernel_Point_Std_Dev_Matrix(ii, jj) = RMSE_VE;
        end
    end
    resultKernel{dtc} = Kernel_Point_Std_Dev_Matrix;


    % test smooth parameters for TPS
    TPS_Point_Std_Dev_Matrix = [];
    for jj = 1 : length(TPSsmoothParams)
        [points_std_dev, RMSE_VE] = runTrainingTest_TPSGPA (DataCell, VisVecCell, TestPointsDataCell, TestPointsVisVecCell, ...
            TPSsmoothParams(jj));
        TPS_Point_Std_Dev_Matrix(jj) = RMSE_VE;
    end    
    resultTPS{dtc} = TPS_Point_Std_Dev_Matrix;


    handles = [];
    figure("Position", [100, 100, 400, 500]); hold on;
    for ii = 1 : length(kernelParams)
        h  = plot(smoothParams, Kernel_Point_Std_Dev_Matrix(ii, :));
        handles = [handles, h];
    end
    legend(handles, kernelParamsStr);    
    set(gca, "XScale", "log");
    hold off;
    title("Kernel GPA")


    figure("Position", [600, 100, 400, 500]); hold on;
    plot(TPSsmoothParams, TPS_Point_Std_Dev_Matrix);
    set(gca, "XScale", "log");
    hold off;
    title("TPS GPA");

    pause(10);
end


%save result_tuning_params.mat resultKernel resultTPS













function [points_std_dev, RMSE_VE] = runTrainingTest_KernelGPA (DataCell, VisVecCell, TestPointsDataCell, TestPointsVisVecCell, ...
    smooth_param, kernel_param)

flagJointlyOptimizeLamdaRigidTransformation = true;

% Kernel-GPA
GPA_KERNEL = KernelGPA;
GPA_KERNEL.verbose = false;
GPA_KERNEL.smoothParam = smooth_param;

%GPA_KERNEL.kernelScaleParam= kernel_param;
GPA_KERNEL.quantilePercentParam = kernel_param;

GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
GPA_KERNEL.run (flagJointlyOptimizeLamdaRigidTransformation);


% transform the point-clouds to the reference frame
for ii = 1 : length(TestPointsDataCell)
    data = GPA_KERNEL.transformPoints(TestPointsDataCell{ii}, ii);
    transTestDataCell{ii} = data;
end

[points_std_dev, mean_template, RMSE_VE] = ComputePointwiseConsistency (transTestDataCell, TestPointsVisVecCell);

end




function [points_std_dev, RMSE_VE] = runTrainingTest_TPSGPA (DataCell, VisVecCell, TestPointsDataCell, TestPointsVisVecCell, ...
    smoothParam)

flagJointlyOptimizeLamdaRigidTransformation = true;

% TPS-GPA. TPS(5)
GPA_TPS = DefGPA('TPS');
GPA_TPS.verbose = false;
GPA_TPS.dimFeatureSpace =5^3;
GPA_TPS.addDataByArray(DataCell, VisVecCell);
GPA_TPS.run (smoothParam, flagJointlyOptimizeLamdaRigidTransformation);
GPA_TPS.rsdError();

% transform the point-clouds to the reference frame
for ii = 1 : length(TestPointsDataCell)
    data = GPA_TPS.transformPoints(TestPointsDataCell{ii}, ii);
    transTestDataCell{ii} = data;
end
                
[points_std_dev, mean_template, RMSE_VE] = ComputePointwiseConsistency (transTestDataCell, TestPointsVisVecCell);

end





