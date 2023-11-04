clear all
close all
clc

addpath('../DefGPA', '../');
addpath('../SGPAv1.0', '../SGPAv1.0/util');

centerePerAxis = [3, 5, 7];

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end



figformat = 'png';

datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug'};
bestParams = {10, 100, 100, 0.01, 0.01};

n = [2,3];

figPos =  [0, 0, 250, 200];
%figPos =  [0, 0, 150, 120];

imgDPI = 300;


% first component move image rightwards
% second component move image downwards
originArray = {[350; 300], [600; 450], [700; 450]};

CropRectArr = {[240, 210, 200, 200],  [70, 50, 1040-1, 800-1],  [165, 50, 1000-1, 880-1]};


for dchoice = [2, 3]

dataName = datasetCell{dchoice};
    
[DataCell, VisVecCell, ImgFileCell] = readDataset (dataName);

smoothParam = bestParams{dchoice};

cropRectangle = CropRectArr{dchoice};


GPA_EUC = GPA_Functor('EUC-ALL');
GPA_EUC.addDataByArray(DataCell, VisVecCell);
GPA_EUC.run();
GlobalRotation =  getGlobalRotation (dataName, GPA_EUC.mShape);
[~, varImg, ~] = IntensityVisualization (GPA_EUC, ImgFileCell, GlobalRotation, originArray{dchoice});
fprintf(1, 'image before crop: %d, %d\n', size(varImg,1), size(varImg,2));
varImg = imcrop(uint8(varImg), cropRectangle);
fprintf(1, 'image after crop: %d, %d\n', size(varImg,1), size(varImg,2));
varImg = imresize(varImg, 0.4);
fprintf(2, 'image resized: %d, %d\n', size(varImg,1), size(varImg,2));
varImg = rescale(varImg);
fig_intensity_std = figure(1);
imshow(varImg);
ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1;
fig_intensity_std.Position = figPos;
fig_saved_name = [fdir, 'Intensity_Deviation_', dataName, '_EUC'];
imwrite(varImg,  [fig_saved_name, '.png']);
exportgraphics(fig_intensity_std, [fig_saved_name, '.pdf'], 'ContentType','Vector');
pause(0.5);




GPA_AFF = GPA_Functor('AFF-ALL');
GPA_AFF.addDataByArray(DataCell, VisVecCell);
GPA_AFF.run();
GlobalRotation =  getGlobalRotation (dataName, GPA_AFF.mShape);
[~, varImg, ~] = IntensityVisualization (GPA_AFF, ImgFileCell, GlobalRotation, originArray{dchoice});
varImg = imcrop(uint8(varImg), cropRectangle);
varImg = imresize(varImg, 0.4);
varImg = rescale(varImg);
fig_intensity_std = figure(2);
imshow(varImg);
ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1;
fig_intensity_std.Position = figPos;
fig_saved_name = [fdir, 'Intensity_Deviation_', dataName, '_AFF'];
imwrite(varImg,  [fig_saved_name, '.png']);
exportgraphics(fig_intensity_std, [fig_saved_name, '.pdf'], 'ContentType','Vector');
pause(0.5);




GPA = DefGPA('AFFINE');
GPA.flagJointlyOptimizeLamdaRigidTransformation = false;
GPA.dimFeatureSpace = 0;
GPA.addDataByArray(DataCell, VisVecCell);
GPA.run (smoothParam);
GlobalRotation =  getGlobalRotation (dataName, GPA.mShape);
[~, varImg, ~]  = IntensityVisualization (GPA, ImgFileCell, GlobalRotation, originArray{dchoice});
varImg = imcrop(uint8(varImg), cropRectangle);
varImg = imresize(varImg, 0.4);
varImg = rescale(varImg);
clear GPA;
fig_intensity_std = figure(3);
imshow(varImg);
ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1;
fig_intensity_std.Position =figPos;
fig_saved_name = [fdir, 'Intensity_Deviation_', dataName, '_AFFINE'];
imwrite(varImg,  [fig_saved_name, '.png']);
exportgraphics(fig_intensity_std, [fig_saved_name, '.pdf'], 'ContentType','Vector');
pause(0.5);


for kk = 1 : size(centerePerAxis,2)

GPA = DefGPA('TPS');
GPA.flagJointlyOptimizeLamdaRigidTransformation = false;
GPA.dimFeatureSpace = centerePerAxis(kk)^2;
GPA.addDataByArray(DataCell, VisVecCell);
GPA.run (smoothParam);
GlobalRotation =  getGlobalRotation (dataName, GPA.mShape);
[~, varImg, ~]  = IntensityVisualization (GPA, ImgFileCell, GlobalRotation, originArray{dchoice});
varImg = imcrop(uint8(varImg), cropRectangle);
varImg = imresize(varImg, 0.4);
varImg = rescale(varImg);
clear GPA;
fig_intensity_std = figure(kk+3);
imshow(varImg);
ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1;
fig_intensity_std.Position = figPos;
fig_saved_name = [fdir, 'Intensity_Deviation_', dataName, '_TPS_', num2str(centerePerAxis(kk))];
imwrite(varImg,  [fig_saved_name, '.png']);
exportgraphics(fig_intensity_std, [fig_saved_name, '.pdf'], 'ContentType','Vector');
pause(0.5);

end

end










function Rot =  getGlobalRotation (dataName, refShape)

[DataCell, VisVecCell] = readDataset (dataName);

GPA_EUC_DAT = GPA_Functor('EUC-ALL');
GPA_EUC_DAT.addDataByArray(DataCell, VisVecCell);
GPA_EUC_DAT.run();

Rot = OrthoProcrustes (refShape, GPA_EUC_DAT.mShape);
end


function [R, t] = OrthoProcrustes (D1, D2)

meanD1 = mean(D1, 2);

meanD2 = mean(D2, 2);

M = (D1 - meanD1) * (D2 - meanD2)';

[U, ~, V] = svd(M);

R = V * U';

if det(R) < 0
    
    fprintf(2, '\n There exists reflection between solutions! \n');
    
end

t = R * meanD1 - meanD2;

end









