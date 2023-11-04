clear all
close all
clc


fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end


figformat = 'pdf';

figSize = [100, 50, 250, 250];


addpath('../DefGPA', '../');
addpath('../SGPAv1.0', '../SGPAv1.0/util');

centerePerAxis = [3, 5, 7];
dimFeatureSpace = {centerePerAxis.^2,  centerePerAxis.^2, centerePerAxis.^2, centerePerAxis.^3, centerePerAxis.^3, 2*centerePerAxis.^2};



datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug'}

bestParams = {10, 100, 100, 0.1, 0.01, 0.01};




n =[1, 2, 3, 4, 6]

% n = [4, 5, 6];
% 
% n = [4, 6];

%n = 4;



for dchoice = n

[DataCell, VisVecCell] = readDataset (datasetCell{dchoice});

smoothParam = bestParams{dchoice};

numPointClouds = max (size(DataCell));


if(1)

% method: EUC-ALL: EUC_d
EUC_ALL = GPA_Functor('EUC-ALL');
EUC_ALL.addDataByArray(DataCell, VisVecCell);
EUC_ALL.run();
EUC_ALL.rsdError();
% for ii = 1 : numPointClouds
%     tformD_EUC_ALL{ii} = EUC_ALL.transformPoints (DataCell{ii}, ii);
% end
if(size(DataCell{1}, 2) < 500)
    tformD_EUC_ALL = LeaveNCrossValidation (DataCell, VisVecCell, -1, 0, 1, EUC_ALL.mShape);
else
    tformD_EUC_ALL = LeaveNCrossValidation (DataCell, VisVecCell, -1, 0, 40, EUC_ALL.mShape);
end

% method: AFF-ALL:  AFF_d
AFF_ALL = GPA_Functor('AFF-ALL');
AFF_ALL.addDataByArray(DataCell, VisVecCell);
AFF_ALL.run();
AFF_ALL.rsdError();
if(size(DataCell{1}, 2) < 500)
    tformD_AFF_ALL = LeaveNCrossValidation (DataCell, VisVecCell, -1, -1, 1, AFF_ALL.mShape);
else
    tformD_AFF_ALL = LeaveNCrossValidation (DataCell, VisVecCell, -1, -1, 40, AFF_ALL.mShape);
end

% AFFINE_r
GPA_AFFINE = DefGPA('AFFINE');
GPA_AFFINE.flagJointlyOptimizeLamdaRigidTransformation = false;
for ii = 1 : size(DataCell,2)
    GPA_AFFINE.addPointCloud (DataCell{ii}, VisVecCell{ii});    
end
GPA_AFFINE.run ();
GPA_AFFINE.rsdError();
if(size(DataCell{1}, 2) < 500)
    tformD_AFFINE = LeaveNCrossValidation (DataCell, VisVecCell, 0, inf, 1, GPA_AFFINE.mShape);
else
    tformD_AFFINE = LeaveNCrossValidation (DataCell, VisVecCell, 0, inf, 40, GPA_AFFINE.mShape);
end


GPA_TPS_Case1 = DefGPA('TPS');
GPA_TPS_Case1.flagJointlyOptimizeLamdaRigidTransformation = false;
for ii = 1 : size(DataCell,2)
    GPA_TPS_Case1.addPointCloud (DataCell{ii}, VisVecCell{ii}); 
end
GPA_TPS_Case1.dimFeatureSpace = dimFeatureSpace{dchoice}(1);
GPA_TPS_Case1.run (smoothParam);
GPA_TPS_Case1.rsdError();
if(size(DataCell{1}, 2) < 500)
    tformD_TPS_Case1 = LeaveNCrossValidation (DataCell, VisVecCell, dimFeatureSpace{dchoice}(1), smoothParam, 1, GPA_TPS_Case1.mShape);
else
    tformD_TPS_Case1 = LeaveNCrossValidation (DataCell, VisVecCell, dimFeatureSpace{dchoice}(1), smoothParam, 40, GPA_TPS_Case1.mShape);
end


GPA_TPS_Case2 = DefGPA('TPS');
GPA_TPS_Case2.flagJointlyOptimizeLamdaRigidTransformation = false;
for ii = 1 : size(DataCell,2)
    GPA_TPS_Case2.addPointCloud (DataCell{ii}, VisVecCell{ii}); 
end
GPA_TPS_Case2.dimFeatureSpace = dimFeatureSpace{dchoice}(2);
GPA_TPS_Case2.run (smoothParam);
GPA_TPS_Case2.rsdError();
if(size(DataCell{1}, 2) < 500)
    tformD_TPS_Case2 = LeaveNCrossValidation (DataCell, VisVecCell, dimFeatureSpace{dchoice}(2), smoothParam, 1, GPA_TPS_Case2.mShape);
else
    tformD_TPS_Case2 = LeaveNCrossValidation (DataCell, VisVecCell, dimFeatureSpace{dchoice}(2), smoothParam, 40, GPA_TPS_Case2.mShape);
end


GPA_TPS_Case3 = DefGPA('TPS');
GPA_TPS_Case3.flagJointlyOptimizeLamdaRigidTransformation = false;
for ii = 1 : size(DataCell,2)
    GPA_TPS_Case3.addPointCloud (DataCell{ii}, VisVecCell{ii}); 
end
GPA_TPS_Case3.dimFeatureSpace = dimFeatureSpace{dchoice}(3);
GPA_TPS_Case3.run (smoothParam);
GPA_TPS_Case3.rsdError();
if(size(DataCell{1}, 2) < 500)
    tformD_TPS_Case3 = LeaveNCrossValidation (DataCell, VisVecCell, dimFeatureSpace{dchoice}(3), smoothParam, 1, GPA_TPS_Case3.mShape);
else
    tformD_TPS_Case3 = LeaveNCrossValidation (DataCell, VisVecCell, dimFeatureSpace{dchoice}(3), smoothParam, 40, GPA_TPS_Case3.mShape);
end

if (dchoice == 5)
    save(['exp_Landmark_visualization_', datasetCell{dchoice}], 'EUC_ALL', 'AFF_ALL', 'GPA_AFFINE', 'GPA_TPS_Case1', 'GPA_TPS_Case2', 'GPA_TPS_Case3', 'tformD_EUC_ALL', 'tformD_AFF_ALL', 'tformD_AFFINE', 'tformD_TPS_Case1', 'tformD_TPS_Case2', 'tformD_TPS_Case3');
end

else

if (dchoice == 5)    
    load(['exp_Landmark_visualization_', datasetCell{dchoice}], 'EUC_ALL', 'AFF_ALL', 'GPA_AFFINE', 'GPA_TPS_Case1', 'GPA_TPS_Case2', 'GPA_TPS_Case3', 'tformD_EUC_ALL', 'tformD_AFF_ALL', 'tformD_AFFINE', 'tformD_TPS_Case1', 'tformD_TPS_Case2', 'tformD_TPS_Case3');
end
    
end


%- - ---- visualization

hFig1 = figure(100 * dchoice + 1);
set(hFig1, 'Position', figSize);
plotLandmarkRegistration (100 * dchoice + 1, EUC_ALL.mShape, tformD_EUC_ALL, VisVecCell, EUC_ALL.mShape);
hFig1 = tightfig(hFig1);
%set(gca,'LooseInset',get(gca,'TightInset'));
fig_saved_name = [fdir,  'exp_Landmark_visualization_', datasetCell{dchoice}, '_EUC_d'];
%saveas (hFig1, fig_saved_name, figformat);
print(hFig1,  '-dpdf', fig_saved_name, '-painters');
%ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat],'Resolution',300);
close all

hFig2 = figure(100 * dchoice + 2);
set(hFig2, 'Position', figSize);
plotLandmarkRegistration (100 * dchoice + 2, AFF_ALL.mShape, tformD_AFF_ALL, VisVecCell, EUC_ALL.mShape);
hFig2 = tightfig(hFig2);
fig_saved_name = [fdir,  'exp_Landmark_visualization_', datasetCell{dchoice}, '_AFF_d'];
%saveas (hFig2, fig_saved_name, figformat);
print(hFig2,  '-dpdf', fig_saved_name, '-painters');
%ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat],'Resolution',300);
close all

hFig3 = figure(100 * dchoice + 3);
set(hFig3, 'Position', figSize);
plotLandmarkRegistration (100 * dchoice + 3, GPA_AFFINE.mShape, tformD_AFFINE, VisVecCell, EUC_ALL.mShape);
hFig3 = tightfig(hFig3);
fig_saved_name = [fdir,  'exp_Landmark_visualization_', datasetCell{dchoice}, '_AFF_r'];
%saveas (hFig3, fig_saved_name, figformat);
print(hFig3,  '-dpdf', fig_saved_name, '-painters');
%ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat],'Resolution',300);
close all

hFig4 = figure(100 * dchoice + 4);
set(hFig4, 'Position', figSize);
plotLandmarkRegistration (100 * dchoice + 4, GPA_TPS_Case1.mShape, tformD_TPS_Case1, VisVecCell, EUC_ALL.mShape)
hFig4 = tightfig(hFig4);
fig_saved_name = [fdir,  'exp_Landmark_visualization_', datasetCell{dchoice}, '_TPS_r_', num2str(centerePerAxis(1))];
%saveas (hFig4, fig_saved_name, figformat);
print(hFig4,  '-dpdf', fig_saved_name, '-painters');
%ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat],'Resolution',300);
close all

hFig5 = figure(100 * dchoice + 5);
set(hFig5, 'Position', figSize);
plotLandmarkRegistration (100 * dchoice + 5, GPA_TPS_Case2.mShape, tformD_TPS_Case2, VisVecCell, EUC_ALL.mShape)
hFig5 = tightfig(hFig5);
fig_saved_name = [fdir,  'exp_Landmark_visualization_', datasetCell{dchoice}, '_TPS_r_', num2str(centerePerAxis(2))];
%saveas (hFig5, fig_saved_name, figformat);
print(hFig5,  '-dpdf', fig_saved_name, '-painters');
%ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat],'Resolution',300);
close all

hFig6 = figure(100 * dchoice + 6);
set(hFig6, 'Position', figSize);
plotLandmarkRegistration (100 * dchoice + 6, GPA_TPS_Case3.mShape, tformD_TPS_Case3, VisVecCell, EUC_ALL.mShape)
hFig6 = tightfig(hFig6);
fig_saved_name = [fdir,  'exp_Landmark_visualization_', datasetCell{dchoice}, '_TPS_r_', num2str(centerePerAxis(3))];
%saveas (hFig6, fig_saved_name, figformat);
print(hFig6,  '-dpdf', fig_saved_name, '-painters');
%ax = gca; exportgraphics(ax, [fig_saved_name, '.', figformat],'Resolution',300);
%close all

end



% The solution to the similarity Procrustes problem
% sR, t = minimize || sR * D1  + t  - D2 ||
function [sR, t] = SimilairtyProcrustes (D1, D2)

meanD1 = mean(D1, 2);

meanD2 = mean(D2, 2);

M = (D1 - meanD1) * (D2 - meanD2)';

N = (D1 - meanD1) * (D1 - meanD1)';

[U, ~, V] = svd(M);

R = V * U';

if det(R) < 0
    
    fprintf(2, '\nThere exists reflection between solutions! \n');
    
end

s = trace(R * M) / trace(N);

sR = s * R;

t =  meanD2 - sR * meanD1;

end



function [R, t] = EuclideanProcrustes (D1, D2)

meanD1 = mean(D1, 2);

meanD2 = mean(D2, 2);

M = (D1 - meanD1) * (D2 - meanD2)';

[U, ~, V] = svd(M);

R = V * U';

if det(R) < 0
    
    fprintf(2, '\n There exists reflection between solutions! \n');
    
end

t =  meanD2 - R * meanD1;

end



function plotLandmarkRegistration (fid, refShape, transDataCell, VisVecCell, gaugeRef)

n = max(size(transDataCell));
m = size(refShape, 2);
dim = size(refShape, 1);

% get cross validation error
normalizer = 0; cve = 0;

for ii = 1 : n
    
    pointsDeviation = (transDataCell{ii} - refShape) .* VisVecCell{ii};
    
    sqrtError =  norm(pointsDeviation, 'fro');
    
    cve = cve + sqrtError * sqrtError;
    
    normalizer = normalizer + nnz(VisVecCell{ii});
    
end

cve = sqrt(cve/(normalizer));


% visualization

% put everything in the coordinate frame of the gaugeRef shape.
% [R, t] = EuclideanProcrustes(refShape, gaugeRef);
[R, t] = SimilairtyProcrustes(refShape, gaugeRef);

refShape = R * refShape + t;
for ii = 1 : n
    transDataCell{ii} = R * transDataCell{ii} + t;
end


% compute position covariance

positionCov = zeros(dim, m);

for jj = 1 : m
    
    PointsN = zeros(dim, n);
    
    for  ii = 1 : n
        
        PointsN(:, ii) = transDataCell{ii}(:, jj);
        
    end
    
    pr = refShape(:, jj);
    
    PointsNError = PointsN - pr;
    
    for dd = 1 : dim
        positionCov(dd, jj) = norm (PointsNError(dd, :), 'fro') / sqrt(n);
    end
    
end


% visualzation:
% 2D point clouds: plot out all landmark positions directly
% 3D point clouds: plot the mean(reference shape), and the position deviations, color encoded.

figure(fid);

% visualization of 2D point clouds

if dim == 2
    
    hold on
    
    DT = delaunayTriangulation(refShape');
    
    for ii = 1 : n
        
        D = transDataCell{ii};
        
        TRI = DT.ConnectivityList;
        
        TRIh = hideNode (TRI, VisVecCell{ii});
        
        triplot (TRIh, D(1,:), D(2,:),  '-b.',  'MarkerSize', 10);

    end
    
    TRI = DT.ConnectivityList;
    
    triplot (TRI, refShape(1,:), refShape(2,:),  '-g.', 'LineWidth', 1.5,  'MarkerSize', 20);
    
    axis image;
    
    hold off
    
    axis equal; grid on;
    
    hf=get(gcf, 'CurrentAxes');
    set(hf, 'YDir', 'Reverse');
    
    title(['CVE = ', num2str(cve, '%.2f')], 'horizontalAlignment', 'center', 'FontSize', 12);
end


% visualization of 3D points clouds

if dim == 3
    codeVal = mean(positionCov, 1);
    axis equal; grid on;
    if (floor(fid/100) == 4)
        scatter3(refShape(1,:),refShape(2,:),refShape(3,:),  10*codeVal, codeVal, 'filled');
        view(-10, 10);
         caxis([0, 20]);
        cb = colorbar('FontSize',10, 'Location', 'northoutside');
        %colormap jet
        colormap turbo
        x1=get(gca,'position');
        pos = get(cb,'Position');
        pos(4) = pos(4)/3;
        set(cb,'position', pos, 'AxisLocation', 'in');
        set(gca,'position',x1);
%         xticks([-50, 0, 50, 100]);
%         zticks([-50, 0, 50]);
        title(['CVE = ', num2str(cve, '%.2f')], 'horizontalAlignment', 'center', 'FontSize', 12);
    end
    if (floor(fid/100) == 5)
        scatter3(refShape(1,:),refShape(2,:),refShape(3,:),  codeVal, codeVal, 'filled');
        caxis([0, 9]);
        colormap turbo
        xlim([-100, 150]);
        ylim([-80, 120]);
        zlim([-120, 100]);
        view(-30, 30);
        cb = colorbar('FontSize',10, 'Location', 'eastoutside');
        x1=get(gca,'position');
        pos = get(cb,'Position');
        pos(3) = pos(3)/3;
        set(cb,'position', pos, 'AxisLocation', 'in');
        set(gca,'position',x1);
        title(['CVE = ', num2str(cve, '%.2f')], 'horizontalAlignment', 'center', 'FontSize', 8);
    end
    if (floor(fid/100) == 6)
        scatter3(refShape(1,:),refShape(2,:),refShape(3,:),  30*codeVal, codeVal, 'filled');
        view(85, 5);
        cb = colorbar('FontSize',10, 'Location', 'northoutside');
        %colormap jet
        colormap turbo
        x1=get(gca,'position');
        pos = get(cb,'Position');
        pos(4) = pos(4)/3;
        set(cb,'position', pos, 'AxisLocation', 'in');
        set(gca,'position',x1);
        title(['CVE = ', num2str(cve, '%.2f')], 'horizontalAlignment', 'center', 'FontSize', 12);
    end
%    xlabel('x'); ylabel('y'); zlabel('z');
end








end






function TRIh = hideNode (TRI, VIS)

TRIh=[];

for i= 1 : size(TRI,1)
    if VIS(TRI(i,1)) && VIS(TRI(i,2)) && VIS(TRI(i,3))
    
        TRIh = [TRIh; TRI(i,:)]; 
    
    end
end

end


