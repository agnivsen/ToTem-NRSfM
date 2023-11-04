
clear all
close all
clc

addpath('../DefGPA', '../');

% ----- Dataset Information ---

datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug', 'TOPACS'};


DefGPA_TPS_bestParams = {10, 100, 100, 0.1, 0.01, 0.01};
dimensions = {2, 2, 2, 3, 3, 3, 3}; 


KGPA_bestQuantilePercentParams = [0.2,  0.2,  0.2,  0.2,  0.2,  0.2];
KGPA_bestSmoothParams = [0.05, 0.05, 0.05, 0.05,  0.05, 0.05];


flagJointlyOptimizeLamdaRigidTransformation = true;

% By setting Nfold to infinity, the method
% run_N_Fold_CrossValidation(Nfold)
% generates Leave-1-out cross-validation.
Nfold = inf;




for dchoice = [1,2,3,4]



datasetName = datasetCell{dchoice};
[DataCell, VisVecCell] = readDataset (datasetName);


GPA_RIGID = RigidGPA;
GPA_RIGID.addDataByArray(DataCell, VisVecCell);
GPA_RIGID.run();
[cve_rigid,trans_pts_rigid] = GPA_RIGID.run_N_Fold_CrossValidation(Nfold);
GPA_RIGID


GPA_AFFINE = DefGPA('AFFINE');
GPA_AFFINE.dimFeatureSpace = 0;
GPA_AFFINE.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_AFFINE.addDataByArray(DataCell, VisVecCell);
GPA_AFFINE.run (DefGPA_TPS_bestParams{dchoice}, flagJointlyOptimizeLamdaRigidTransformation);
GPA_AFFINE.rsdError();
[cve_affine,trans_pts_affine] = GPA_AFFINE.run_N_Fold_CrossValidation(Nfold);
GPA_AFFINE


GPA_TPS3 = DefGPA('TPS');
GPA_TPS3.dimFeatureSpace =3^dimensions{dchoice};
GPA_TPS3.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_TPS3.addDataByArray(DataCell, VisVecCell);
GPA_TPS3.run (DefGPA_TPS_bestParams{dchoice}, flagJointlyOptimizeLamdaRigidTransformation);
GPA_TPS3.rsdError();
[cve_tps3,trans_pts_tps3] = GPA_TPS3.run_N_Fold_CrossValidation(Nfold);
GPA_TPS3


GPA_TPS5 = DefGPA('TPS');
GPA_TPS5.dimFeatureSpace =5^dimensions{dchoice};
GPA_TPS5.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_TPS5.addDataByArray(DataCell, VisVecCell);
GPA_TPS5.run (DefGPA_TPS_bestParams{dchoice}, flagJointlyOptimizeLamdaRigidTransformation);
GPA_TPS5.rsdError();
[cve_tps5,trans_pts_tps5] = GPA_TPS5.run_N_Fold_CrossValidation(Nfold);
GPA_TPS5


GPA_TPS7 = DefGPA('TPS');
GPA_TPS7.dimFeatureSpace = 7^dimensions{dchoice};
GPA_TPS7.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_TPS7.addDataByArray(DataCell, VisVecCell);
GPA_TPS7.run (DefGPA_TPS_bestParams{dchoice}, flagJointlyOptimizeLamdaRigidTransformation);
GPA_TPS7.rsdError();
[cve_tps7,trans_pts_tps7] = GPA_TPS7.run_N_Fold_CrossValidation(Nfold);
GPA_TPS7


GPA_KERNEL = KernelGPA;
GPA_KERNEL.smoothParam = KGPA_bestSmoothParams(dchoice);
GPA_KERNEL.quantilePercentParam= KGPA_bestQuantilePercentParams(dchoice);
GPA_KERNEL.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
GPA_KERNEL.run ();
GPA_KERNEL.rsdError();
[cve_kernel,trans_pts_kernel] = GPA_KERNEL.run_N_Fold_CrossValidation(Nfold);
GPA_KERNEL




%%

gaugeShape = GPA_RIGID.mShape;

% rotate point-clouds to the same reference frame
rect_R = OrhoProcurestes(GPA_RIGID.mShape, gaugeShape);
for ii = 1 : length(trans_pts_rigid)
    trans_pts_rigid{ii} = rect_R * trans_pts_rigid{ii};
end

rect_R = OrhoProcurestes(GPA_AFFINE.mShape, gaugeShape);
for ii = 1 : length(trans_pts_affine)
    trans_pts_affine{ii} = rect_R * trans_pts_affine{ii};
end

rect_R = OrhoProcurestes(GPA_TPS3.mShape, gaugeShape);
for ii = 1 : length(trans_pts_tps3)
    trans_pts_tps3{ii} = rect_R * trans_pts_tps3{ii};
end

rect_R = OrhoProcurestes(GPA_TPS5.mShape, gaugeShape);
for ii = 1 : length(trans_pts_tps5)
    trans_pts_tps5{ii} = rect_R * trans_pts_tps5{ii};
end

rect_R = OrhoProcurestes(GPA_TPS7.mShape, gaugeShape);
for ii = 1 : length(trans_pts_tps7)
    trans_pts_tps7{ii} = rect_R * trans_pts_tps7{ii};
end

rect_R = OrhoProcurestes(GPA_KERNEL.mShape, gaugeShape);
for ii = 1 : length(trans_pts_kernel)
    trans_pts_kernel{ii} = rect_R * trans_pts_kernel{ii};
end

%%



fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end


close all;
% We use correspondences that occur in all 6 point-clouds to validate registeration results
figSize = [100, 50, 400, 400];


tfontsize = 16;

tfig1 = figure(1);
stdvar1 = VisualizePointCloudsVariation (trans_pts_rigid, VisVecCell);
tfig1.Name = ['Rigid (', num2str(stdvar1, '%.2f'), ')'];
if (dchoice == 4)
view(-100, 10);
end
axis equal tight; pause(0.1); axis off;
set(tfig1, 'Position', figSize);
title(['std-dev= ', num2str(stdvar1, '%.2f')], 'Interpreter', 'Latex', 'FontSize', tfontsize);
exportgraphics(tfig1, [fdir, datasetName, '_GPA_Visualization_Rigid', '.pdf'], 'ContentType','vector');


tfig2 = figure(2);
stdvar2 = VisualizePointCloudsVariation (trans_pts_affine, VisVecCell);
tfig2.Name = ['Affine (', num2str(stdvar2, '%.2f'), ')'];
if (dchoice == 4)
view(-100, 10);
end
axis equal tight; pause(0.1); axis off;
set(tfig2, 'Position', figSize);
title(['std-dev = ', num2str(stdvar2, '%.2f')], 'Interpreter', 'Latex', 'FontSize', tfontsize);
exportgraphics(tfig2, [fdir, datasetName, '_GPA_Visualization_Affine', '.pdf'], 'ContentType','vector');


tfig3 = figure(3);
stdvar3 = VisualizePointCloudsVariation (trans_pts_tps3, VisVecCell);
tfig3.Name = ['TPS(3) (', num2str(stdvar3, '%.2f'), ')'];
if (dchoice == 4)
view(-100, 10); axis equal tight;
end
axis equal tight; pause(0.1); axis off;
set(tfig3, 'Position', figSize);
title(['std-dev = ', num2str(stdvar3, '%.2f')], 'Interpreter', 'Latex', 'FontSize', tfontsize);
exportgraphics(tfig3, [fdir, datasetName, '_GPA_Visualization_TPS3', '.pdf'], 'ContentType','vector');


tfig4 = figure(4);
stdvar4 = VisualizePointCloudsVariation (trans_pts_tps5, VisVecCell);
tfig4.Name = ['TPS(5) (', num2str(stdvar4, '%.2f'), ')'];
if (dchoice == 4)
view(-100, 10); axis equal tight;
end
axis equal tight; pause(0.1); axis off;
set(tfig4, 'Position', figSize);
title(['std-dev = ', num2str(stdvar4, '%.2f')], 'Interpreter', 'Latex', 'FontSize', tfontsize);
exportgraphics(tfig4, [fdir, datasetName, '_GPA_Visualization_TPS5', '.pdf'], 'ContentType','vector');


tfig5 = figure(5);
stdvar5 = VisualizePointCloudsVariation (trans_pts_tps7, VisVecCell);
tfig5.Name = ['TPS(7) (', num2str(stdvar5, '%.2f'), ')'];
if (dchoice == 4)
view(-100, 10);
end
axis equal tight; pause(0.1); axis off;
set(tfig5, 'Position', figSize);
title(['std-dev = ', num2str(stdvar5, '%.2f')], 'Interpreter', 'Latex', 'FontSize', tfontsize);
exportgraphics(tfig5, [fdir, datasetName, '_GPA_Visualization_TPS7', '.pdf'], 'ContentType','vector');


tfig6 = figure(6);
stdvar6 = VisualizePointCloudsVariation (trans_pts_kernel, VisVecCell);
tfig6.Name = ['Kernel (', num2str(stdvar6, '%.2f'), ')'];
if (dchoice == 4)
view(-100, 10);
end
axis equal tight; pause(0.1); axis off;
set(tfig6, 'Position', figSize);
title(['std-dev = ', num2str(stdvar6, '%.2f')], 'Interpreter', 'Latex', 'FontSize', tfontsize);
exportgraphics(tfig6, [fdir, datasetName, '_GPA_Visualization_Kernel', '.pdf'], 'ContentType','vector');


fprintf(1, 'GPA - Rigid = %f\n', stdvar1);
fprintf(1, 'GPA - Affine = %f\n', stdvar2);
fprintf(1, 'GPA - TPS(5) = %f\n', stdvar3);
fprintf(1, 'GPA - TPS(7) = %f\n', stdvar4);
fprintf(1, 'GPA - Kernel = %f\n', stdvar6);


end





function signs = GetRectSigns (Shape1, ShapeGauge)

v = sum(Shape1, 2);
vg = sum(ShapeGauge, 2);

ptc = 2;

v =  Shape1(: ,ptc);
vg = ShapeGauge(: ,ptc);

signs = sign(v) ./ sign(vg);


end


% orthogonal Procrustes for  R * D1 = D2
function R = OrhoProcurestes (D1, D2)

meanD1 = mean(D1, 2);

meanD2 = mean(D2, 2);

M = (D1 - meanD1) * (D2 - meanD2)';

[U, ~, V] = svd(M);

R = V * U';

end



