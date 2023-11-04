clear all;
close all;
clc;

%-%
addpath('../DefGPA', '../');

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end



% red = [0.55    0    0.];
% redMat = [0.6353    0.0784    0.1843];
% 
% blue = [0    0    0.55];
% lightblue = [0,  0.4,  0.95];
% 
% green = [0    0.55    0];
% brightgreeen = [0.3137    0.9412    0.0863];
% 
% shinyred = [0.65, 0, 0];
% shinygreen = [0, 0.65, 0];
% shinyblue = [0, 0, 0.65];
%  
% yellow = [0.9, 0.7, 0];
% yellowMat = [0.9294    0.6941    0.1255];
% 
% grassyellow = [0.6314    0.6314         0];
% 
% purple = [0.4941    0.1843    0.5569];
% 
% brightcycan = [0    0.7    1];
% darkcyan = [0,  0.5843,  0.7137];
% 
% brightorange = [1.0000    0.4118    0.1608];
% 
% brown = [0.6000    0.4510    0.1020];
% 
% gray = [0.7, 0.7, 0.7];
% 
% deepmegenta = [0.7608    0.1294    0.5176];
% 
% darkmegenta = [0.6902    0.2000    0.5020];



colors = {
[ 0.8500    0.3250    0.0980 ], ...
[ 0.9290    0.6940    0.1250 ], ...
[ 0.4940    0.1840    0.5560 ], ...
[ 0         0         0      ], ...
}


colors = {
[ 0    0    0 ], ...
[ 0    0    0 ], ...
[ 0    0    0 ], ...
[ 0    0    0 ], ...
}


colors = {
[ 0    0.8    0 ], ...
[ 0    0.8    0 ], ...
[ 0    0.8    0 ], ...
[ 0    0.8    0 ], ...
}



colors = {
[ 0.9    0.7    0 ], ...
[ 0.9    0.7    0 ], ...
[ 0.9    0.7    0 ], ...
[ 0.9    0.7    0 ], ...
}






length = 4;
linewidth = 2;
pos = [0.2, 0.2, 0.2] + length;
fontsize = 25;


fdir = '../../PaperDraft/figure/';


load('xyzPoints');
ptCloud = pointCloud(xyzPoints);

%ptCloud = pcread('teapot.ply');
pts = ptCloud.Location';

pts = pts + [4; 4; 0];

num_points = size(pts, 2);


%C = turbo (num_points);
C = winter (num_points);
C = flipud(winter(num_points));

C = uint8(round(C*255));

C1 = uint8(ones(num_points, 1) * [0, 0, 0]); 
C2 = uint8(ones(num_points, 1) * [0, 0, 255]);
C3 = uint8(ones(num_points, 1) * [255, 0, 0]);



hFig1 = figure;
pc=pointCloud(pts');
pc.Color = C1;
pcshow(pc)
xlabel('x');ylabel('y');zlabel('z');
hold on;
ax1 = quiver3(0, 0, 0, length, 0, 0, 'LineWidth', linewidth, 'Color', colors{1});
ax2 = quiver3(0, 0, 0, 0, length, 0, 'LineWidth', linewidth, 'Color', colors{2});
ax3 = quiver3(0, 0, 0, 0, 0, length, 'LineWidth', linewidth, 'Color', colors{3});
text(pos(1), 0, 0, 'X', 'FontSize', fontsize, 'Color', ax1.Color);
text(0, pos(2), 0, 'Y', 'FontSize', fontsize, 'Color', ax2.Color);
text(0, 0, pos(3), 'Z', 'FontSize', fontsize, 'Color', ax3.Color);
text(0.0, 0.5, -0.5, 'O', 'FontSize', fontsize, 'Color', colors{4});
hold off;
axis equal off;
exportgraphics(gca, [fdir, 'example1.pdf'], 'ContentType','vector');



center = mean(pts, 2);
centered_pts = pts - center;



centered_pc=pointCloud(centered_pts');
centered_pc.Color = C2;
hFig2 = figure;
pcshow(centered_pc);
xlabel('x');ylabel('y');zlabel('z');
hold on;
ax1 = quiver3(0, 0, 0, length, 0, 0, 'LineWidth', linewidth, 'Color', colors{1});
ax2 = quiver3(0, 0, 0, 0, length, 0, 'LineWidth', linewidth, 'Color', colors{2});
ax3 = quiver3(0, 0, 0, 0, 0, length, 'LineWidth', linewidth, 'Color', colors{3});
text(pos(1), 0, 0, 'X', 'FontSize', fontsize, 'Color', ax1.Color);
text(0, pos(2), 0, 'Y', 'FontSize', fontsize, 'Color', ax2.Color);
text(0, 0, pos(3), 'Z', 'FontSize', fontsize, 'Color', ax3.Color);
text(0.0, 0.5, -0.5, 'O', 'FontSize', fontsize, 'Color', colors{4});
hold off;
axis equal off;
exportgraphics(gca, [fdir, 'example2.pdf'], 'ContentType','vector');



pc_cov = double(centered_pts * centered_pts.');
[Q, S] = eig(pc_cov);
R = Q';
rotated_centered_pts = R * centered_pts;


rotated_centered_pc=pointCloud(rotated_centered_pts');
rotated_centered_pc.Color = C3;
hFig3 = figure;
pcshow(rotated_centered_pc);
xlabel('x');ylabel('y');zlabel('z');
hold on;
ax1 = quiver3(0, 0, 0, length, 0, 0, 'LineWidth', linewidth, 'Color', colors{1});
ax2 = quiver3(0, 0, 0, 0, length, 0, 'LineWidth', linewidth, 'Color', colors{2});
ax3 = quiver3(0, 0, 0, 0, 0, length, 'LineWidth', linewidth, 'Color', colors{3});
text(pos(1), 0, 0, 'X', 'FontSize', fontsize, 'Color', ax1.Color);
text(0, pos(2), 0, 'Y', 'FontSize', fontsize, 'Color', ax2.Color);
text(0, 0, pos(3), 'Z', 'FontSize', fontsize, 'Color', ax3.Color);
text(0.0, 0.5, -0.5, 'O', 'FontSize', fontsize, 'Color', colors{4});
hold off;
axis equal off;
exportgraphics(gca, [fdir, 'example3.pdf'], 'ContentType','vector');



final_test = S - rotated_centered_pts * rotated_centered_pts'



function plotAxis (color_info, headstyle, linewidth)

color_info = 'r';
headWidth = 4;
headLength = 4;
linewidth = 1;
headstyle = 'vback2';

ah = annotation('arrow',...
    'headStyle',headstyle,'HeadLength',headLength,'HeadWidth',headWidth, 'Color', color_info, "LineWidth", linewidth);
set(ah,'parent',gca);
%set(ah, 'Position', [t(1), t(2), R(1, 1), R(2, 1)]);



ah = annotation('arrow',...
    'headStyle',headstyle,'HeadLength',headLength,'HeadWidth',headWidth, 'Color', color_info , "LineWidth", linewidth);
set(ah,'parent',gca);
set(ah, 'Position', [t(1), t(2), R(1, 2), R(2, 2)]);

end
