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


figformat = 'pdf';

datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug'};
bestParams = {10, 100, 100, 0.1, 0.1};

n = 1;

figPos =  [100, 50, 200, 200];

cropRectangle = [240, 210, 200-1, 200-1];

cropOrigin = [cropRectangle(1); cropRectangle(2)];


referMarker = 'go';   color1 = 'g';
tformMarker = 'bd';  color2 = 'b';




n_chosen_img = 5;

markersize = 4;


% first component move image rightwards
% second component move image downwards
originArray = {[350; 300], [600; 450], [700; 450]};

for dchoice = n

dataName = datasetCell{dchoice};
    
[DataCell, VisVecCell, ImgFileCell] = readDataset (dataName);

smoothParam = bestParams{dchoice};




GPA_EUC = GPA_Functor('EUC-ALL');
GPA_EUC.addDataByArray(DataCell, VisVecCell);
GPA_EUC.run();
GlobalRotation =  getGlobalRotation (dataName, GPA_EUC.mShape);
[tformImgCell, ~, ~] = IntensityVisualization (GPA_EUC, ImgFileCell, GlobalRotation, originArray{dchoice});
tformLandmarks = GlobalRotation * GPA_EUC.transformPoints (DataCell{n_chosen_img}, n_chosen_img) + [originArray{dchoice}(1); originArray{dchoice}(2)] - cropOrigin;
refShape = GlobalRotation * GPA_EUC.mShape + [originArray{dchoice}(1); originArray{dchoice}(2)] - cropOrigin;
fig_tform_imgs = figure;
cropImg = imcrop(tformImgCell{n_chosen_img}, cropRectangle);
imshow(uint8(cropImg));
hold on
plot(refShape(1, :), refShape(2, :), referMarker, 'MarkerSize', 1.0 * markersize, 'markerfacecolor',color1);
plot(tformLandmarks(1, logical(VisVecCell{n_chosen_img})), tformLandmarks(2, logical(VisVecCell{n_chosen_img})), tformMarker, 'MarkerSize', markersize, 'markerfacecolor',color2);
hold off
ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1;
fig_tform_imgs.Position = figPos;
fig_tform_imgs = tightfig(fig_tform_imgs);
fig_saved_name = [fdir, 'transformed_imgs_landmarks_', dataName, '_EUC'];
%saveas (fig_tform_imgs, fig_saved_name, figformat );
%print(fig_tform_imgs,  '-dpdf', fig_saved_name, '-painters', '-loose');
ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat], 'ContentType','Vector');
pause(0.5);
clear GPA_EUC;




GPA_AFF = GPA_Functor('AFF-ALL');
GPA_AFF.addDataByArray(DataCell, VisVecCell);
GPA_AFF.run();
GlobalRotation =  getGlobalRotation (dataName, GPA_AFF.mShape);
[tformImgCell, ~, ~] = IntensityVisualization (GPA_AFF, ImgFileCell, GlobalRotation, originArray{dchoice});
tformLandmarks = GlobalRotation * GPA_AFF.transformPoints (DataCell{n_chosen_img}, n_chosen_img) + [originArray{dchoice}(1); originArray{dchoice}(2)] - cropOrigin;
refShape = GlobalRotation * GPA_AFF.mShape + [originArray{dchoice}(1); originArray{dchoice}(2)] - cropOrigin;
fig_tform_imgs = figure;
cropImg = imcrop(tformImgCell{n_chosen_img}, cropRectangle);
imshow(uint8(cropImg));
hold on
plot(refShape(1, :), refShape(2, :), referMarker, 'MarkerSize', 1.0 * markersize, 'markerfacecolor', color1);
plot(tformLandmarks(1, logical(VisVecCell{n_chosen_img})), tformLandmarks(2, logical(VisVecCell{n_chosen_img})), tformMarker, 'MarkerSize', markersize, 'markerfacecolor',color2);
hold off
ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1;
fig_tform_imgs.Position = figPos;
fig_tform_imgs = tightfig(fig_tform_imgs);
fig_saved_name = [fdir, 'transformed_imgs_landmarks_', dataName, '_AFF'];
%saveas (fig_tform_imgs, fig_saved_name, figformat );
%print(fig_tform_imgs,  '-dpdf', fig_saved_name, '-painters');
ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat], 'ContentType','Vector');
pause(0.5);
clear GPA_AFF;





GPA = DefGPA('AFFINE');
GPA.flagJointlyOptimizeLamdaRigidTransformation = false;
GPA.dimFeatureSpace = 0;
GPA.addDataByArray(DataCell, VisVecCell);
GPA.run (smoothParam);
GlobalRotation =  getGlobalRotation (dataName, GPA.mShape);
[tformImgCell, ~, ~]  = IntensityVisualization (GPA, ImgFileCell, GlobalRotation, originArray{dchoice});
tformLandmarks = GlobalRotation * GPA.transformPoints (DataCell{n_chosen_img}, n_chosen_img) + [originArray{dchoice}(1); originArray{dchoice}(2)] - cropOrigin;
refShape = GlobalRotation * GPA.mShape + [originArray{dchoice}(1); originArray{dchoice}(2)] - cropOrigin;
fig_tform_imgs = figure;
cropImg = imcrop(tformImgCell{n_chosen_img}, cropRectangle);
imshow(uint8(cropImg));
hold on
plot(refShape(1, :), refShape(2, :), referMarker, 'MarkerSize', 1.0 * markersize, 'markerfacecolor', color1);
plot(tformLandmarks(1, logical(VisVecCell{n_chosen_img})), tformLandmarks(2, logical(VisVecCell{n_chosen_img})), tformMarker, 'MarkerSize', markersize, 'markerfacecolor',color2);
hold off
ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1;
fig_tform_imgs.Position =figPos;
fig_tform_imgs = tightfig(fig_tform_imgs);
fig_saved_name = [fdir, 'transformed_imgs_landmarks_', dataName, '_AFFINE'];
%saveas (fig_tform_imgs, fig_saved_name, figformat );
%print(fig_tform_imgs,  '-dpdf', fig_saved_name, '-painters');
ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat], 'ContentType','Vector');
pause(0.5);
clear GPA;


for kk = 1 : size(centerePerAxis,2)

GPA = DefGPA('TPS');
GPA.flagJointlyOptimizeLamdaRigidTransformation = false;
GPA.dimFeatureSpace = centerePerAxis(kk)^2;
GPA.addDataByArray(DataCell, VisVecCell);
GPA.run (smoothParam);
GlobalRotation =  getGlobalRotation (dataName, GPA.mShape);
[tformImgCell, ~, ~]  = IntensityVisualization (GPA, ImgFileCell, GlobalRotation, originArray{dchoice});
tformLandmarks = GlobalRotation * GPA.transformPoints (DataCell{n_chosen_img}, n_chosen_img) + [originArray{dchoice}(1); originArray{dchoice}(2)] - cropOrigin;
refShape = GlobalRotation * GPA.mShape + [originArray{dchoice}(1); originArray{dchoice}(2)] - cropOrigin;
fig_tform_imgs = figure;
cropImg = imcrop(tformImgCell{n_chosen_img}, cropRectangle);
imshow(uint8(cropImg));
hold on
plot(refShape(1, :), refShape(2, :), referMarker, 'MarkerSize', 1.0 * markersize, 'markerfacecolor',color1);
plot(tformLandmarks(1, logical(VisVecCell{n_chosen_img})), tformLandmarks(2, logical(VisVecCell{n_chosen_img})), tformMarker, 'MarkerSize', markersize, 'markerfacecolor',color2);
hold off
ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1;
fig_tform_imgs.Position = figPos;
fig_tform_imgs = tightfig(fig_tform_imgs);
fig_saved_name = [fdir, 'transformed_imgs_landmarks_', dataName, '_TPS_', num2str(centerePerAxis(kk))];
%saveas (fig_tform_imgs, fig_saved_name, figformat );
%print(fig_tform_imgs,  '-dpdf', fig_saved_name, '-painters');
ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat], 'ContentType','Vector');
pause(0.5);
clear GPA;

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









