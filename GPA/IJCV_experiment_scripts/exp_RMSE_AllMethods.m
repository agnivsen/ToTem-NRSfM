clear all
close all
clc

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end



addpath('../DefGPA', '../');
addpath('../SGPAv1.0', '../SGPAv1.0/util');

centerePerAxis = [3, 5, 7];



            
fid = fopen ([fdir, 'exp_landmark_residual_table.tex'], 'w');
fprintf (fid, '\\begin{tabular}{@{} llcccccc @{}}\n');
fprintf (fid, '\\toprule \n');
fprintf (fid, 'Dataset &   &  $\\ast$EUC\\_d & $\\ast$AFF\\_d & AFF\\_r & TPS\\_r($%d^d$) & TPS\\_r($%d^d$) & TPS\\_r($%d^d$)  \\\\ \n', centerePerAxis(1), centerePerAxis(2), centerePerAxis(3));
fprintf (fid, '\\midrule \n');

formatSpec_RMSE_d = ['%s ', ' & RMSE\\_d & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\ \n'];
formatSpec_RMSE_r = [' ', ' & RMSE\\_r & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\ \n'];
formatSpec_CVE   = [' ', ' & CVE & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\ \n'];
formatSpec_Time   = [' ', ' & Time & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n'];



datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug'};
datasetNameCell = {'Face', 'Bag', 'Pillow', 'LiTS', 'Liver', 'ToyRug'};
dimFeatureSpace = {centerePerAxis.^2,  centerePerAxis.^2, centerePerAxis.^2, centerePerAxis.^3, centerePerAxis.^3, 2*centerePerAxis.^2};

bestParams = {10, 100, 100, 0.1, 0.01, 0.01};


for ii = 1 : size(datasetCell,2)

[DataCell, VisVecCell] = readDataset (datasetCell{ii});

[EUC_d, AFF_d, AFF_r, TPS_r_case1, TPS_r_case2, TPS_r_case3] = runAllMethods(DataCell, VisVecCell, bestParams{ii}, dimFeatureSpace{ii});
RMSE_d = [EUC_d.RMSE_d, AFF_d.RMSE_d, AFF_r.RMSE_d, TPS_r_case1.RMSE_d, TPS_r_case2.RMSE_d, TPS_r_case3.RMSE_d];
RMSE_r = [EUC_d.RMSE_r, AFF_d.RMSE_r, AFF_r.RMSE_r, TPS_r_case1.RMSE_r, TPS_r_case2.RMSE_r, TPS_r_case3.RMSE_r];
CVE = [EUC_d.CVE, AFF_d.CVE, AFF_r.CVE, TPS_r_case1.CVE, TPS_r_case2.CVE, TPS_r_case3.CVE];
Time = [EUC_d.time, AFF_d.time, AFF_r.time, TPS_r_case1.time, TPS_r_case2.time, TPS_r_case3.time];

fprintf (fid, formatSpec_RMSE_d, datasetNameCell{ii}, RMSE_d);
fprintf (fid, formatSpec_RMSE_r, RMSE_r);
fprintf (fid, formatSpec_CVE, CVE);
fprintf (fid, formatSpec_Time, Time);

end


fprintf (fid, '\\bottomrule \n');
fprintf (fid, '\\end{tabular}');
fclose (fid);






function [EUC_d, AFF_d, AFF_r, TPS_r_case1, TPS_r_case2, TPS_r_case3] = runAllMethods(DataCell, VisVecCell, smoothParam, dimFeatureSpace)

% method: EUC-ALL
EUC_ALL = GPA_Functor('EUC-ALL');
EUC_ALL.addDataByArray(DataCell, VisVecCell);
EUC_ALL.run();
EUC_ALL.rsdError();

EUC_d.RMSE_r = EUC_ALL.rsd_ref;
EUC_d.RMSE_d = EUC_ALL.rsd_dat;
EUC_d.time = EUC_ALL.time(1);
EUC_d.CVE  = computeCrossValidationError(DataCell, VisVecCell, -1, 0, EUC_ALL.mShape);
clear EUC_ALL


% method: AFF-ALL
AFF_ALL = GPA_Functor('AFF-ALL');
AFF_ALL.addDataByArray(DataCell, VisVecCell);
AFF_ALL.run();
AFF_ALL.rsdError();

AFF_d.RMSE_r = AFF_ALL.rsd_ref;
AFF_d.RMSE_d = AFF_ALL.rsd_dat;
AFF_d.time = AFF_ALL.time(1);
AFF_d.CVE  = computeCrossValidationError(DataCell, VisVecCell, -1, -1, AFF_ALL.mShape);
clear AFF_ALL


% AFFINE
GPA_AFFINE = DefGPA('AFFINE');
GPA_AFFINE.flagJointlyOptimizeLamdaRigidTransformation = false;
GPA_AFFINE.addDataByArray(DataCell, VisVecCell);
GPA_AFFINE.run ();
GPA_AFFINE.rsdError();

AFF_r.RMSE_r = GPA_AFFINE.rsd_ref;
AFF_r.RMSE_d = GPA_AFFINE.rsd_dat;
AFF_r.time = GPA_AFFINE.time(1);
AFF_r.CVE  = computeCrossValidationError(DataCell, VisVecCell, 0, inf, GPA_AFFINE.mShape);

clear GPA_AFFINE



% TPS warp, centerePerAxis(1)^d
GPA_TPS = DefGPA('TPS');
GPA_TPS.flagJointlyOptimizeLamdaRigidTransformation = false;
GPA_TPS.addDataByArray(DataCell, VisVecCell);
GPA_TPS.dimFeatureSpace = dimFeatureSpace(1);
GPA_TPS.run (smoothParam);
GPA_TPS.rsdError();

TPS_r_case1.RMSE_r = GPA_TPS.rsd_ref;
TPS_r_case1.RMSE_d = GPA_TPS.rsd_dat;
TPS_r_case1.time = GPA_TPS.time(1);
TPS_r_case1.CVE = computeCrossValidationError(DataCell, VisVecCell, dimFeatureSpace(1), smoothParam, GPA_TPS.mShape);

clear GPA_TPS



% TPS warp, centerePerAxis(2)^d
GPA_TPS = DefGPA('TPS');
GPA_TPS.flagJointlyOptimizeLamdaRigidTransformation = false;
GPA_TPS.addDataByArray(DataCell, VisVecCell);
GPA_TPS.dimFeatureSpace = dimFeatureSpace(2);
GPA_TPS.run (smoothParam);
GPA_TPS.rsdError();

TPS_r_case2.RMSE_r = GPA_TPS.rsd_ref;
TPS_r_case2.RMSE_d = GPA_TPS.rsd_dat;
TPS_r_case2.time = GPA_TPS.time(1);
TPS_r_case2.CVE = computeCrossValidationError(DataCell, VisVecCell, dimFeatureSpace(2), smoothParam, GPA_TPS.mShape);

clear GPA_TPS




% TPS warp, centerePerAxis(3)^d
GPA_TPS = DefGPA('TPS');
GPA_TPS.flagJointlyOptimizeLamdaRigidTransformation = false;
GPA_TPS.addDataByArray(DataCell, VisVecCell);
GPA_TPS.dimFeatureSpace = dimFeatureSpace(3);
GPA_TPS.run (smoothParam);
GPA_TPS.rsdError();

TPS_r_case3.RMSE_r = GPA_TPS.rsd_ref;
TPS_r_case3.RMSE_d = GPA_TPS.rsd_dat;
TPS_r_case3.time = GPA_TPS.time(1);
TPS_r_case3.CVE = computeCrossValidationError(DataCell, VisVecCell, dimFeatureSpace(3), smoothParam, GPA_TPS.mShape);

clear GPA_TPS



end