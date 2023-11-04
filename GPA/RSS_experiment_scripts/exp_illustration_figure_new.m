clear all;
close all;
clc;

addpath('../DefGPA', '../');

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end


rng(0);

% simulate ground-truth landmarks
dxy = 1;
gt_x = [1 2 3 1 1 2 3 1 1 1];
gt_y = [6 6 6 5 4 4 4 3 2 1];
gt_points = dxy * [gt_x;  gt_y];    gt_center = dxy* [2; 4];
gt_points = gt_points - gt_center;
% simulate ground-truth robot poses
drbt = 6;
pos1 = drbt * [-1; 0];
pos2 = drbt * [0; -1];
pos3 = drbt * [1;  0];
pos4 = drbt * [0;  1];
rot1 = [1, 1;  -1, 1] ./ sqrt(2);
rot2 = [1, -1;  1, 1] ./ sqrt(2);
rot3 = [-1, -1; 1, -1] ./ sqrt(2);
rot4 = [-1, 1; -1, -1] ./ sqrt(2);
% simulate visibility
vis1 = true(1, 10);  vis1([2, 3, 7]) = false;
vis2 = true(1, 10);  vis2([1, 2, 3]) = false;
vis3 = true(1, 10);  vis3([1, 4, 9, 10]) = false;
vis4 = true(1, 10);  vis4([8, 9, 10]) = false;

sigma_deform = 0;
points1 = gt_points +  sigma_deform * SimulateSmoothDeformationField (gt_points, 10, 100);
points2 = gt_points +  sigma_deform * SimulateSmoothDeformationField (gt_points, 10, 100);
points3 = gt_points +  sigma_deform * SimulateSmoothDeformationField (gt_points, 10, 100);
points4 = gt_points +  sigma_deform * SimulateSmoothDeformationField (gt_points, 10, 100);

sigma_noise = 0.1;
% relative measurements
% Global coordinate system  [i, j] at [0, 0]
% Local robot coordinate system  [e1, e2] = [i, j] * [r1, r2] at [c1, c2]
local_points1 = rot1' * (points1 - pos1) + sigma_noise * min(max(randn(2, length(points1)), -3), 3);
local_points2 = rot2' * (points2 - pos2) + sigma_noise * min(max(randn(2, length(points2)), -3), 3);
local_points3 = rot3' * (points3 - pos3) + sigma_noise * min(max(randn(2, length(points3)), -3), 3);
local_points4 = rot4' * (points4 - pos4) + sigma_noise * min(max(randn(2, length(points4)), -3), 3);


% Create Data into Cell
DataCell{1} = local_points1; VisVecCell{1} = vis1;
DataCell{2} = local_points2; VisVecCell{2} = vis2;
DataCell{3} = local_points3; VisVecCell{3} = vis3;
DataCell{4} = local_points4; VisVecCell{4} = vis4;


%%

GPA_RIGID = RigidGPA;
GPA_RIGID.addDataByArray(DataCell, VisVecCell);
GPA_RIGID.run();

GPA_KERNEL = KernelGPA;
GPA_KERNEL.smoothParam = 0.1;
GPA_KERNEL.quantilePercentParam= 0.05;
GPA_KERNEL.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
GPA_KERNEL.run ();



%%


fontsize = 7;
axismargin = 0.5;
markersize = 60;

gt_headstyle = 'cback1';
est_headstyle = 'vback2';


marker0 = 'o';
marker5 = 'o';

marker1 = '>';
marker2 = '^';
marker3 = '<';
marker4 = 'v';

% 
% marker1 = 's';
% marker2 = '+';
% marker3 = '*';
% marker4 = 'x';
 
% color1 = '#0E9DFF'
% color2 = '#FF0000'
% color3 = '#800080'
% color4 = '#FFA500'
 
% color1 = '#B0E188'
% color2 = '#2077B5'
% color3 = '#05B9C7'
% color4 = '#A8CBE4'
 
% color1 = '#313BD0'
% color2 = '#9A22F8'
% color3 = '#00F2F2'
% color4 = '#00B2FC'




color0 = '#979797'
color1 =  'r'
color2 = '#7F1EAC'
color3 = '#21A576'
color4 = '#F8A427'
color5 ='k'


color1 = '#3742FA';
color2 = '#05C46B';
color3 = [0 1 1];
color4 = '#FC5C65';
color3 = '#18C1FF';


axiscolor0 = color0;
axiscolor1 = color1;
axiscolor2 = color2;
axiscolor3 = color3;
axiscolor4 = color4;
axiscolor5 = color5;


titlecolor0 = axiscolor0;
titlecolor1 = 'k';
titlecolor2 = 'k';
titlecolor3 = 'k';
titlecolor4 = 'k';
titlecolor5 = 'k';



gt_posecolor1 = axiscolor1;
gt_posecolor2 = axiscolor2;
gt_posecolor3 = axiscolor3;
gt_posecolor4 = axiscolor4;

est_posecolor1 = axiscolor1;
est_posecolor2 = axiscolor2;
est_posecolor3 = axiscolor3;
est_posecolor4 = axiscolor4;


hFig = figure("Name", "Point-Clouds", "Position", [0, 200, 700, 250]);
tfig = tiledlayout(4, 10, "TileSpacing","tight");


axlinewidth = 0.5;
axfontsize = 6;
linewidth = 1;



ax0 = nexttile(1, [4, 4]);
ax0.XColor = axiscolor0; ax0.YColor = axiscolor0;
ax0.FontSize = axfontsize;
ax0.LineWidth = axlinewidth;
hold on;
scatter(gt_points(1, :), gt_points(2, :), markersize*0.5, 'filled', 'Marker', marker0, "MarkerEdgeColor", color0, "MarkerFaceColor", color0);
plotPose([rot1, pos1], 'pose1', gt_posecolor1, gt_headstyle, 1.2);
plotPose([rot2, pos2], 'pose2', gt_posecolor2, gt_headstyle, 1.2);
plotPose([rot3, pos3], 'pose3', gt_posecolor3, gt_headstyle, 1.2);
plotPose([rot4, pos4], 'pose4', gt_posecolor4, gt_headstyle, 1.2);
hold off;
axis equal; box on;
axlimx = ax0.XLim; axlimy = ax0.YLim;
ax0.XLim = [axlimx(1)-axismargin,   axlimx(2)+axismargin];
ax0.YLim = [axlimy(1)-axismargin,   axlimy(2)+axismargin];
%title("ground-truth point-cloud and poses", 'FontSize', fontsize, 'FontWeight','normal', 'Color', titlecolor0);




% transformed point-clouds
ax1 = nexttile(5);
ax1.XColor = axiscolor1; ax1.YColor = axiscolor1;
ax1.FontSize = axfontsize;
ax1.LineWidth = axlinewidth;
hold on;
hpt1 = scatter(local_points1(1, vis1),local_points1(2, vis1), markersize*0.5, 'Marker', marker1, "MarkerEdgeColor", color1, "MarkerFaceColor", 'none', 'LineWidth', linewidth);
hold off;
axis equal square;
box on;
%title("point-cloud 1", 'FontSize', fontsize, 'FontWeight','normal', 'Color', titlecolor1);

ax2 = nexttile(15);
ax2.XColor = axiscolor2; ax2.YColor = axiscolor2;
ax2.FontSize = axfontsize;
ax2.LineWidth = axlinewidth;
hold on;
hpt2 = scatter(local_points2(1, vis2), local_points2(2, vis2), markersize*0.5, 'Marker', marker2, "MarkerEdgeColor", color2, "MarkerFaceColor", 'none', 'LineWidth', linewidth);
hold off;
axis equal square; box on;
%title("point-cloud 2", 'FontSize', fontsize, 'FontWeight','normal', 'Color', titlecolor2);

ax3 = nexttile(25);
ax3.XColor = axiscolor3; ax3.YColor = axiscolor3;
ax3.FontSize = axfontsize;
ax3.LineWidth = axlinewidth;
hold on;
hpt3 = scatter(local_points3(1, vis3), local_points3(2, vis3), markersize*0.5, 'Marker', marker3, "MarkerEdgeColor", color3, "MarkerFaceColor", 'none', 'LineWidth', linewidth);
hold off;
axis equal square; box on;
%title("point-cloud 3", 'FontSize', fontsize, 'FontWeight','normal', 'Color', titlecolor3);

ax4 = nexttile(35);
ax4.XColor = axiscolor4; ax4.YColor = axiscolor4;
ax4.FontSize = axfontsize;
ax4.LineWidth = axlinewidth;
hold on;
hpt4 = scatter(local_points4(1, vis4), local_points4(2, vis4), markersize*0.5, 'Marker', marker4, "MarkerEdgeColor", color4, "MarkerFaceColor", 'none', 'LineWidth', linewidth);
hold off;
axis equal square; box on;
%title("point-cloud 4", 'FontSize', fontsize, 'FontWeight','normal', 'Color', titlecolor4);





% Registeration result Rigid

%GPA_HANDLE = GPA_RIGID;
GPA_HANDLE = GPA_KERNEL;
refpointcloud = GPA_HANDLE.mShape;
pose_estimates = GPA_HANDLE.mPoses;
transform_points1 = GPA_HANDLE.transformPoints(DataCell{1}, 1);
transform_points2 = GPA_HANDLE.transformPoints(DataCell{2}, 2);
transform_points3 = GPA_HANDLE.transformPoints(DataCell{3}, 3);
transform_points4 = GPA_HANDLE.transformPoints(DataCell{4}, 4);

ax5 = nexttile(7,  [4, 4]);
ax5.XColor = axiscolor5; ax5.YColor = axiscolor5;
ax5.FontSize = axfontsize;
ax5.LineWidth = axlinewidth;
hold on;
hr1 = scatter(transform_points1(1, vis1), transform_points1(2, vis1), markersize, 'Marker', marker1, "MarkerEdgeColor", color1, "MarkerFaceColor", 'none', 'LineWidth', linewidth);
hr2 = scatter(transform_points2(1, vis2), transform_points2(2, vis2), markersize, 'Marker', marker2, "MarkerEdgeColor", color2, "MarkerFaceColor", 'none', 'LineWidth', linewidth);
hr3 = scatter(transform_points3(1, vis3), transform_points3(2, vis3), markersize, 'Marker', marker3, "MarkerEdgeColor", color3, "MarkerFaceColor", 'none', 'LineWidth', linewidth);
hr4 = scatter(transform_points4(1, vis4), transform_points4(2, vis4), markersize, 'Marker', marker4, "MarkerEdgeColor", color4, "MarkerFaceColor", 'none', 'LineWidth', linewidth);
hrr = scatter(refpointcloud(1, :), refpointcloud(2, :), markersize*0.1, 'Marker', marker5, "MarkerEdgeColor", color5, "MarkerFaceColor", color5);
plotPose(pose_estimates{1}, 'pose1', est_posecolor1, est_headstyle);
plotPose(pose_estimates{2}, 'pose2', est_posecolor2, est_headstyle);
plotPose(pose_estimates{3}, 'pose3', est_posecolor3, est_headstyle);
plotPose(pose_estimates{4}, 'pose4', est_posecolor4, est_headstyle);
hold off;
axis equal; box on;
axlimx = ax5.XLim; axlimy = ax5.YLim;
ax5.XLim = [axlimx(1)-axismargin,   axlimx(2)+axismargin];
ax5.YLim = [axlimy(1)-axismargin,   axlimy(2)+axismargin];
%title({"global registration by GPA"}, 'FontSize', fontsize, 'FontWeight','normal', 'Color', titlecolor5);


lgd = legend([hpt1, hpt2, hpt3, hpt4], {'point-cloud1', 'point-cloud2', 'point-cloud3', 'point-cloud4'}, 'box', 'off', 'FontSize', 8, 'Orientation','horizontal');
lgd.Layout.Tile = 'south';

%lgd = legend([hr1, hr2, hr3, hr4, hrr], {'point-cloud1', 'point-cloud2', 'point-cloud3', 'point-cloud4', 'estimated point-cloud'}, 'box', 'off', ...
 %   'Orientation','vertical', 'FontSize', fontsize);
%lgd.Layout.Tile = "south";
%lgd.Layout.TileSpan= [1, 1];






exportgraphics(hFig, [fdir, 'GPA_Illustration', '.pdf'], 'ContentType','vector');




function plotPose(pose, text_info, color_info, headstyle, linewidth)
if ~exist("info_text", 'var')
    text_info = 'pose';
end
if ~exist("color_info", 'var')
    color_info = 'k';
end
if ~exist("headstyle", "var")
    headstyle = 'cback2';
end
if ~exist("linewidth", "var")
    linewidth = 0.5;
end
t = pose(:, 3);
R = pose(:, [1,2]);
scatter(t(1), t(2), "Marker", '.', 'MarkerEdgeColor',color_info, 'MarkerFaceColor', color_info);

headWidth = 4;
headLength = 4;

ah = annotation('arrow',...
    'headStyle',headstyle,'HeadLength',headLength,'HeadWidth',headWidth, 'Color', color_info, "LineWidth", linewidth);
set(ah,'parent',gca);
set(ah, 'Position', [t(1), t(2), R(1, 1), R(2, 1)]);

% ah = annotation('textarrow', 'string', info_text, 'fontsize', 6, ...
%     'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
ah = annotation('arrow',...
    'headStyle',headstyle,'HeadLength',headLength,'HeadWidth',headWidth, 'Color', color_info , "LineWidth", linewidth);
set(ah,'parent',gca);
set(ah, 'Position', [t(1), t(2), R(1, 2), R(2, 2)]);
end



function smooth_deformation_field =  SimulateSmoothDeformationField (QueryPoints, n_cut_off, nSampleSize)
if ~exist("QueryPoints", 'var') || (exist("QueryPoints", 'var') && isempty(QueryPoints))
    dim = 3;
else
    dim = size(QueryPoints, 1);
end
if ~exist('n_cnt', 'var')
    n_cut_off = 4;
end
if ~exist('nSampleSize', 'var')
    nSampleSize = 100;
end

% Random numbers by Gaussian Distribution
%X = randn(dim, nSampleSize); X = max(X, -3); X = min(X, 3);
% Random numbers by Uniform Distribution
X = 2*(rand(dim, nSampleSize) - 0.5);
% Fourier Transform
F = fft(X, [], 2);
% Cut-off_frequency to filter out high frequency variations
F ( :,  (1+n_cut_off) : (end-n_cut_off+1) ) = 0;
% Inverse Fourier Transform to obtain Smooth Random Numbers
SmoothX = ifft (F, [], 2, 'symmetric');

if exist('QueryPoints', 'var')
    % USE INTERPOLATION FOR QUERY POINTS
    RangeMinMax = [min(QueryPoints, [], 2),  max(QueryPoints, [], 2)];
    lengthXYZ = RangeMinMax (:, 2) - RangeMinMax(:, 1);
    scale = max (lengthXYZ) / (nSampleSize - 1);
    % It is important to USE THE SAME SCALE to ensure uniform deformation in X-Y-Z coordinates
    int_vec = 0 : (nSampleSize-1);
    int_vec(1) = int_vec(1) - 0.1;
    int_vec(end) = int_vec(end) + 0.1;
    vals = scale * ones(dim ,1) * int_vec + RangeMinMax(:, 1);
    smooth_deformation_field = zeros( size(QueryPoints, 1), size(QueryPoints, 2) );

    % obtain the deformation at query-points by interpolation
    for dd = 1 : size(QueryPoints, 1)
        smooth_deformation_field(dd, :) = interp1(vals(dd, :), SmoothX(dd, :), QueryPoints(dd, :));
    end
else
    smooth_deformation_field = SmoothX;
end

if nnz (isnan (smooth_deformation_field))
    fprintf(2, "find error\n");
    smooth_deformation_field = [];
end

end



