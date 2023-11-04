clear all
close all
clc


addpath('../DefGPA', '../');
addpath('../SGPAv1.0', '../SGPAv1.0/util');


sampleChoice = 32

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end


figformat = 'pdf';

dataName = 'stereoToyRug';


centerePerAxis = [3, 5, 7];


figsize = [0, 0, 600, 250];

cropRectangle = [250, 180, 200, 150];




[DataCell, VisVecCell] = readDataset ('StereoToyRug');

smoothParam = 0.01;


EUC_ALL = GPA_Functor('EUC-ALL');
EUC_ALL.addDataByArray(DataCell, VisVecCell);
EUC_ALL.run();
GeneratedShape = EUC_ALL.generatePointCloud (sampleChoice);
hFig1 = visualizeGenenerativelModel (1, sampleChoice, GeneratedShape, cropRectangle);
hFig1.Position = figsize;
hFig1 = tightfig(hFig1);
fig_saved_name = [fdir, 'generated_shape_projection_LR_', dataName, '_EUC'];
%print(hFig1,  '-dpdf', fig_saved_name, '-painters');
exportgraphics(hFig1, [fig_saved_name, '.', figformat],'Resolution',300, 'ContentType','vector');
clear EUC_ALL;


AFF_ALL = GPA_Functor('AFF-ALL');
AFF_ALL.addDataByArray(DataCell, VisVecCell);
AFF_ALL.run();
GeneratedShape = AFF_ALL.generatePointCloud (sampleChoice);
hFig2 = visualizeGenenerativelModel (2, sampleChoice, GeneratedShape, cropRectangle);
hFig2.Position = figsize;
hFig2 = tightfig(hFig2);
fig_saved_name = [fdir, 'generated_shape_projection_LR_', dataName, '_AFF'];
%print(hFig2,  '-dpdf', fig_saved_name, '-painters');
exportgraphics(hFig2, [fig_saved_name, '.', figformat],'Resolution',300, 'ContentType','vector');
clear AFF_ALL


GPA = DefGPA('AFFINE');
GPA.flagJointlyOptimizeLamdaRigidTransformation = false;
GPA.dimFeatureSpace = 0;
GPA.addDataByArray(DataCell, VisVecCell);
GPA.run (smoothParam);
GeneratedShape = GPA.generatePointCloud (sampleChoice);
hFig3 = visualizeGenenerativelModel (3, sampleChoice, GeneratedShape, cropRectangle);
hFig3.Position = figsize;
hFig3 = tightfig(hFig3);
fig_saved_name = [fdir, 'generated_shape_projection_LR_', dataName, '_AFFINE'];
%print(hFig3,  '-dpdf', fig_saved_name, '-painters');
exportgraphics(hFig3, [fig_saved_name, '.', figformat],'Resolution',300, 'ContentType','vector');
clear GPA;




for kk = 1 : size(centerePerAxis,2)

close all;
    
GPA = DefGPA('TPS');
GPA.flagJointlyOptimizeLamdaRigidTransformation = false;
GPA.addDataByArray(DataCell, VisVecCell);
GPA.dimFeatureSpace = 2*centerePerAxis(kk)^2;
GPA.run (smoothParam);
GeneratedShape = GPA.generatePointCloud (sampleChoice);
hFigTPS = visualizeGenenerativelModel (3+kk, sampleChoice, GeneratedShape, cropRectangle);
hFigTPS.Position = figsize;
hFigTPS = tightfig(hFigTPS);
fig_saved_name = [fdir, 'generated_shape_projection_LR_', dataName, '_TPS_', num2str(centerePerAxis(kk))];
%print(hFigTPS,  '-dpdf', fig_saved_name, '-painters');
exportgraphics(hFigTPS, [fig_saved_name, '.', figformat],'Resolution',300, 'ContentType','vector');
clear GPA;

end






function hFig = visualizeGenenerativelModel (fid, sampleChoice, GeneratedShape, cropRectangle)

cropOrigin = [cropRectangle(1); cropRectangle(2)];

allSamples = 1 : 3 : 600;


fdir = '../../dataset/stereoToyRug/';

% number of frames
n = 651;

% number of points
m = 30;

% read left image point tracks
Xl = load([fdir, 'X_left.txt']);

% read right image point tracks
Xr = load([fdir, 'X_right.txt']);

% read 3D points
X = load([fdir, 'X.txt']);

% read camera parameters and reassemble them to projection matrix
Kl = load([fdir, 'Kl.txt']);
Kr = load([fdir, 'Kr.txt']);
Rl = load([fdir, 'Rl.txt']);
Rr = load([fdir, 'Rr.txt']);
Tl = load([fdir, 'Tl.txt']);
Tr = load([fdir, 'Tr.txt']);
Pl = Kl*[Rl  100*Tl];
Pr = Kr*[Rr 100*Tr];


% visualise data sequentially
hFig = figure(fid); clf;
title('red = original points ; blue = reprojections');


for i = allSamples (sampleChoice)
    
    fprintf('frame %03d / %03d\n',i,n);
    
    % read left image
    Il = imread(sprintf([fdir, 'I_left_%03d.jpg'],i));
    
    % read right image
    Ir = imread(sprintf([fdir, 'I_right_%03d.jpg'],i));

    cropIl = imcrop(Il, cropRectangle);

    cropIr = imcrop(Ir, cropRectangle);
    
    % extract left points
    ql = Xl(2*i-1:2*i,:) - cropOrigin;

    % extract right points
    qr = Xr(2*i-1:2*i,:) - cropOrigin;

    % extract 3D points
    Q = 100 * X(3*i-2:3*i,:);

    % reproject 3D points
    hql = Pl*[Q ; ones(1,m)];
    hql = hql(1:2,:)./repmat(hql(3,:),2,1) - cropOrigin;
    hqr = Pr*[Q ; ones(1,m)];
    hqr = hqr(1:2,:)./repmat(hqr(3,:),2,1) - cropOrigin;
    
    

    % reproject 3D points of generated shape
    ghql = Pl*[GeneratedShape ; ones(1,m)];
    ghql = ghql(1:2,:)./repmat(ghql(3,:),2,1) - cropOrigin;
    ghqr = Pr*[GeneratedShape ; ones(1,m)];
    ghqr = ghqr(1:2,:)./repmat(ghqr(3,:),2,1) - cropOrigin;
    

    RMSE_L = norm(ql - ghql, 'fro');
    RMSE_R = norm(qr - ghqr, 'fro');
    
    
    subplot('Position', [0.04, 0.05, 0.45, 0.9]);
    hold off, imshow(cropIl), hold on
    plot(ql(1,:),ql(2,:),'ro','markerfacecolor','r');
    plot(hql(1,:),hql(2,:),'bx');
    plot(ghql(1,:),ghql(2,:),'g*');
    ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1; ax.Title.FontSize = 12;
    title(sprintf('RMSE = %.2f', RMSE_L));
    

    subplot('Position', [0.04+0.02+0.45, 0.05, 0.45, 0.9]);
    hold off, imshow(cropIr), hold on
    plot(qr(1,:),qr(2,:),'ro','markerfacecolor','r');
    plot(hqr(1,:),hqr(2,:),'bx');
    plot(ghqr(1,:),ghqr(2,:),'g*');
    ax = gca;  ax.XAxis.FontSize = 0.1; ax.YAxis.FontSize = 0.1; ax.Title.FontSize = 12;
    title(sprintf('RMSE = %.2f', RMSE_R));
    

    pause(0.1);
    
end

end
