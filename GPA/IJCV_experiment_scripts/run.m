clear all
close all
clc

addpath('../DefGPA', '../');
addpath('../SGPAv1.0', '../SGPAv1.0/util');

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end


exp_Dataset_Sample_Visualization

exp_TPS_Params

exp_Landmark_Visualization

exp_RMSE_AllMethods

exp_Intensity_Visualization
 
exp_Face_visualization

exp_StereoToyRug_Visualziation

exp_Liver_Visualization