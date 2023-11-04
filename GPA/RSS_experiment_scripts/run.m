close all;
clear all;
clc;

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end


addpath('../DefGPA/',  '../');

% this shows kernel parameters for KernelGPA

exp_Kernel_Params_NfoldCVE



exp_TPS_params_NfoldCVE



% this generates error statistic table
exp_RMSE_AllMethods_Kernel_RSS


exp_Lambda_by_MissingPoints


exp_PoseEstimation


exp_TOPACS_visualization


exp_KernelGPA_example_visualization

