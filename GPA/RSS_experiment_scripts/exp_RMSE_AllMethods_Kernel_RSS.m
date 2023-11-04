clear all
close all
clc


fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end


addpath('../DefGPA', '../');

centerePerAxis = [3, 5, 7];


fid = fopen ([fdir, 'exp_landmark_residual_table.tex'], 'w');
fprintf (fid, '\\begin{tabular}{@{} ll|c|cc|cc|cc|cc|cc @{}}\n');
fprintf (fid, '\\toprule \n');
fprintf (fid, '    &   &  EUC & \\multicolumn{2}{c}{AFF} & \\multicolumn{2}{|c}{TPS($%d$)} & \\multicolumn{2}{|c}{TPS($%d$)} & \\multicolumn{2}{|c}{TPS($%d$)}  &  \\multicolumn{2}{|c}{KernelGPA}  \\\\ \n', centerePerAxis(1), centerePerAxis(2), centerePerAxis(3));
fprintf (fid, '\\midrule \n');
fprintf (fid, '\\multicolumn{2}{c|}{$\\boldsymbol{\\Lambda}$ Estimation}   &    &  by \\cite{bai2022ijcv} & proposed & by \\cite{bai2022ijcv} & proposed & by \\cite{bai2022ijcv} & proposed & by \\cite{bai2022ijcv} & proposed  & by \\cite{bai2022ijcv} & proposed  \\\\ \n');


formatSpec_RMSE_d = ['%s ', ' & RMSE\\_d & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f  &  - &  -   \\\\ \n'];
formatSpec_RMSE_r =  ['%s ',  ' & RMSE\\_r & %.2f & %.2f & %.2f & %.2f & %.2f  & %.2f & %.2f & %.2f & %.2f  & %.2f & %.2f  \\\\ \n'];
formatSpec_CVE_Nfold   =      ['%s ',  ' & CVE-$20$fold & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f  & %.2f & %.2f \\\\ \n'];
formatSpec_CVE_Leave1   =      ['%s ',  ' & CVE-Leave$1$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f  & %.2f & %.2f \\\\ \n'];
formatSpec_Time   =    ['%s ', ' & Time [s] & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f  & %.4f & %.4f  \\\\ \n'];


datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug'};
datasetNameCell = {'Face', 'HandBag', 'Pillow', 'LiTS', 'Liver', 'ToyRug'};
dimFeatureSpace = {centerePerAxis.^2,  centerePerAxis.^2, centerePerAxis.^2, centerePerAxis.^3, centerePerAxis.^3, 2*centerePerAxis.^2};


bestParams = {10, 100, 100, 0.1, 0.01, 0.01};



% KGPA_bestSmoothParams = [0.1, 0.1, 0.01, 0.01,  0.003, 0.01];
% KGPA_bestKernelParams = [20, 50,  200,    20, 30,  10];

% new kernel parameters decided by quantiles

KGPA_bestQuantilePercentParams = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
KGPA_bestSmoothParams = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05];


for ii = [1, 2, 3, 4, 6]

[DataCell, VisVecCell] = readDataset (datasetCell{ii});

datasetName = datasetCell{ii};
%load(['training_statistics_params_', datasetName], 'optSmoothParam', 'optKernelParam');

optSmoothParam = KGPA_bestSmoothParams(ii);
optQuantleParam = KGPA_bestQuantilePercentParams(ii);


flagJointlyOptimizeLamdaRigidTransformation = false;
[EUC, AFF_r, TPS_r_case1, TPS_r_case2, TPS_r_case3] = runAllMethods(DataCell, VisVecCell, bestParams{ii}, dimFeatureSpace{ii}, flagJointlyOptimizeLamdaRigidTransformation);
kernel = runKernelGPA(DataCell, VisVecCell, optSmoothParam, optQuantleParam, flagJointlyOptimizeLamdaRigidTransformation);

flagJointlyOptimizeLamdaRigidTransformation = true;
[EUC, AFF_r_new, TPS_r_case1_new, TPS_r_case2_new, TPS_r_case3_new] = runAllMethods(DataCell, VisVecCell, bestParams{ii}, dimFeatureSpace{ii}, flagJointlyOptimizeLamdaRigidTransformation);
kernel_new = runKernelGPA(DataCell, VisVecCell, optSmoothParam, optQuantleParam, flagJointlyOptimizeLamdaRigidTransformation);


RMSE_d = [EUC.RMSE_d, ...
    AFF_r.RMSE_d, AFF_r_new.RMSE_d, ...
    TPS_r_case1.RMSE_d, TPS_r_case1_new.RMSE_d, ...
    TPS_r_case2.RMSE_d, TPS_r_case2_new.RMSE_d, ...
    TPS_r_case3.RMSE_d, TPS_r_case3_new.RMSE_d];
RMSE_r = [EUC.RMSE_r, ...
    AFF_r.RMSE_r, AFF_r_new.RMSE_r, ...
    TPS_r_case1.RMSE_r, TPS_r_case1_new.RMSE_r, ...
    TPS_r_case2.RMSE_r, TPS_r_case2_new.RMSE_r, ...
    TPS_r_case3.RMSE_r, TPS_r_case3_new.RMSE_r, ...
    kernel.RMSE_r, kernel_new.RMSE_r];
CVE_Leave1 = [EUC.CVE_Leave1, ...
    AFF_r.CVE_Leave1, AFF_r_new.CVE_Leave1, ...
    TPS_r_case1.CVE_Leave1, TPS_r_case1_new.CVE_Leave1, ...
    TPS_r_case2.CVE_Leave1, TPS_r_case2_new.CVE_Leave1, ...
    TPS_r_case3.CVE_Leave1, TPS_r_case3_new.CVE_Leave1 ,...
    kernel.CVE_Leave1,  kernel_new.CVE_Leave1];
CVE_Nfold = [EUC.CVE_Nfold, ...
    AFF_r.CVE_Nfold, AFF_r_new.CVE_Nfold, ...
    TPS_r_case1.CVE_Nfold, TPS_r_case1_new.CVE_Nfold, ...
    TPS_r_case2.CVE_Nfold, TPS_r_case2_new.CVE_Nfold, ...
    TPS_r_case3.CVE_Nfold, TPS_r_case3_new.CVE_Nfold,...
    kernel.CVE_Nfold,  kernel_new.CVE_Nfold];
Time = [EUC.time, ...
    AFF_r.time, AFF_r_new.time, ...
    TPS_r_case1.time, TPS_r_case1_new.time, ...
    TPS_r_case2.time, TPS_r_case2_new.time, ...
    TPS_r_case3.time, TPS_r_case3_new.time, ...
    kernel.time, kernel_new.time];

fprintf (fid, '\\midrule \n');
fprintf (fid, formatSpec_RMSE_r, datasetNameCell{ii}, RMSE_r);
fprintf (fid, formatSpec_RMSE_d, ' ', RMSE_d);
fprintf (fid, formatSpec_CVE_Leave1, ' ', CVE_Leave1);
fprintf (fid, formatSpec_CVE_Nfold, ' ', CVE_Nfold);
fprintf (fid, formatSpec_Time, ' ', Time);

end


fprintf (fid, '\\bottomrule \n');
fprintf (fid, '\\end{tabular}');
fclose (fid);




function kernel = runKernelGPA(DataCell, VisVecCell, smoothParam, quantileParam, flagJointlyOptimizeLamdaRigidTransformation)

N = 20;

GPA_KERNEL = KernelGPA;
GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
GPA_KERNEL.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
GPA_KERNEL.smoothParam = smoothParam;
GPA_KERNEL.quantilePercentParam = quantileParam;
GPA_KERNEL.run();

kernel.RMSE_r = GPA_KERNEL.rsdError();
kernel.time = GPA_KERNEL.mTime(1);

kernel.CVE_Nfold = GPA_KERNEL.run_N_Fold_CrossValidation(N);
kernel.CVE_Leave1 = GPA_KERNEL.run_N_Fold_CrossValidation(inf);

end






function [EUC,  AFF_r, TPS_r_case1, TPS_r_case2, TPS_r_case3] = runAllMethods(DataCell, VisVecCell, smoothParam, dimFeatureSpace, flagJointlyOptimizeLamdaRigidTransformation)

N = 20;

% method: EUCLIDEAN
GPA_RIGID = RigidGPA;
GPA_RIGID.addDataByArray(DataCell, VisVecCell);
GPA_RIGID.run();
rsd_r_d = GPA_RIGID.rsdError();

EUC.RMSE_r = rsd_r_d;
EUC.RMSE_d = rsd_r_d;
EUC.time = GPA_RIGID.time(1);

EUC.CVE_Nfold = GPA_RIGID.run_N_Fold_CrossValidation(N);
EUC.CVE_Leave1  = GPA_RIGID.run_N_Fold_CrossValidation(inf);

clear GPA_RIGID



% AFFINE
GPA_AFFINE = DefGPA('AFFINE');
GPA_AFFINE.addDataByArray(DataCell, VisVecCell);
GPA_AFFINE.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
GPA_AFFINE.smoothParam = smoothParam;
GPA_AFFINE.run ();
GPA_AFFINE.rsdError();

AFF_r.RMSE_r = GPA_AFFINE.rsd_ref;
AFF_r.RMSE_d = GPA_AFFINE.rsd_dat;
AFF_r.time = GPA_AFFINE.time(1);

AFF_r.CVE_Nfold = GPA_AFFINE.run_N_Fold_CrossValidation(N);
AFF_r.CVE_Leave1  = GPA_AFFINE.run_N_Fold_CrossValidation(inf);

clear GPA_AFFINE



% TPS warp, centerePerAxis(1)^d
GPA_TPS_CASE1 = DefGPA('TPS');
GPA_TPS_CASE1.addDataByArray(DataCell, VisVecCell);
GPA_TPS_CASE1.dimFeatureSpace = dimFeatureSpace(1);
GPA_TPS_CASE1.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
GPA_TPS_CASE1.smoothParam = smoothParam;
GPA_TPS_CASE1.run ();
GPA_TPS_CASE1.rsdError();

TPS_r_case1.RMSE_r = GPA_TPS_CASE1.rsd_ref;
TPS_r_case1.RMSE_d = GPA_TPS_CASE1.rsd_dat;
TPS_r_case1.time = GPA_TPS_CASE1.time(1);

TPS_r_case1.CVE_Nfold = GPA_TPS_CASE1.run_N_Fold_CrossValidation(N);
TPS_r_case1.CVE_Leave1 = GPA_TPS_CASE1.run_N_Fold_CrossValidation(inf);

clear GPA_TPS_CASE1



% TPS warp, centerePerAxis(2)^d
GPA_TPS_CASE2 = DefGPA('TPS');
GPA_TPS_CASE2.addDataByArray(DataCell, VisVecCell);
GPA_TPS_CASE2.dimFeatureSpace = dimFeatureSpace(2);
GPA_TPS_CASE2.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
GPA_TPS_CASE2.smoothParam = smoothParam;
GPA_TPS_CASE2.run ();
GPA_TPS_CASE2.rsdError();

TPS_r_case2.RMSE_r = GPA_TPS_CASE2.rsd_ref;
TPS_r_case2.RMSE_d = GPA_TPS_CASE2.rsd_dat;
TPS_r_case2.time = GPA_TPS_CASE2.time(1);

TPS_r_case2.CVE_Nfold = GPA_TPS_CASE2.run_N_Fold_CrossValidation(N);
TPS_r_case2.CVE_Leave1 = GPA_TPS_CASE2.run_N_Fold_CrossValidation(inf);

clear GPA_TPS_CASE2



% TPS warp, centerePerAxis(3)^d
GPA_TPS_CASE3 = DefGPA('TPS');
GPA_TPS_CASE3.addDataByArray(DataCell, VisVecCell);
GPA_TPS_CASE3.dimFeatureSpace = dimFeatureSpace(3);
GPA_TPS_CASE3.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
GPA_TPS_CASE3.smoothParam = smoothParam;
GPA_TPS_CASE3.run ();
GPA_TPS_CASE3.rsdError();

TPS_r_case3.RMSE_r = GPA_TPS_CASE3.rsd_ref;
TPS_r_case3.RMSE_d = GPA_TPS_CASE3.rsd_dat;
TPS_r_case3.time = GPA_TPS_CASE3.time(1);

TPS_r_case3.CVE_Nfold = GPA_TPS_CASE3.run_N_Fold_CrossValidation(N);
TPS_r_case3.CVE_Leave1 = GPA_TPS_CASE3.run_N_Fold_CrossValidation(inf);

clear GPA_TPS_CASE3

end