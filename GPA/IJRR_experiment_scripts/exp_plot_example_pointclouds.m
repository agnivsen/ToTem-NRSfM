clear all;
close all;
clc;

%-%
addpath('../DefGPA', '../');

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end



colors = {
[ 0.8500    0.3250    0.0980 ], ...
[ 0.9290    0.6940    0.1250 ], ...
[ 0.4940    0.1840    0.5560 ], ...
}



length = 4;
linewidth = 2;
pos = [0.2, 0.2, 0.2] + length;
fontsize = 25;


fdir = '../../PaperDraft/figure/';


if (0)
%filepath = '../../dataset/bunny/';
%pc = pcread([filepath, 'reconstruction/', 'bun_zipper.ply']);
%pc = pcread([filepath, 'reconstruction/', 'bun_zipper_res3.ply']); % 2,3,4
%pc = pcread([filepath, 'data/', 'bun000.ply']);


%filepath = '../../dataset/dragon/';
%pc = pcread([filepath, 'dragon_recon/', 'dragon_vrip.ply']);

end


filepath = '../../dataset/deformed_2022_11_12/';


colors = {
[ 1    0    0 ], ...
[ 0    0    1 ], ...
[ 0    0.8    0 ], ...
}




ncz = 10;

% cannonical template
xyzpoints_cannonical = readmatrix([filepath, 'karim_liver_0.off'] ,"FileType","text", 'Range', '1:2002');

% deformed point cloud
tmp = pcread([filepath, 'deformed_', num2str(ncz),'.ply']); 
xyzpoints_deformed = tmp.Location;

% deformation field
deform_filed = xyzpoints_deformed - xyzpoints_cannonical;


% observation of visible points
tmp = pcread([filepath, 'visible_', num2str(ncz), '.ply']); 
xyzpoints_observation = tmp.Location;

% visible points
vis = readmatrix([filepath, 'index_', num2str(ncz), '.txt']);
vis_true = false(2002, 1);
vis_true(vis) = true;



cmatrix = ones(size(xyzpoints_cannonical)).* colors{1};
pc_cannonical = pointCloud(xyzpoints_cannonical, 'Color', cmatrix);


xyzpoints_deformed = xyzpoints_deformed(vis_true, :);
cmatrix = ones(size(xyzpoints_deformed)).* colors{2};
pc_deformed = pointCloud(xyzpoints_deformed, 'Color', cmatrix);


cmatrix = ones(size(xyzpoints_observation)).* colors{3};
pc_observation = pointCloud(xyzpoints_observation, 'Color', cmatrix);



figure
hold on;
pcshow(pc_cannonical);
pcshow(pc_deformed);
hold off;


figure
hold on;
pcshow(pc_cannonical);
X = xyzpoints_cannonical(vis_true, 1);
Y = xyzpoints_cannonical(vis_true, 2);
Z = xyzpoints_cannonical(vis_true, 3);
U = deform_filed(vis_true, 1);
V = deform_filed(vis_true, 2);
W = deform_filed(vis_true, 3);
quiver3(X,Y,Z,U,V,W);
pcshow(pc_deformed);
hold off;


figure
hold on;
pcshow(pc_observation)
hold off;


axis equal tight;




