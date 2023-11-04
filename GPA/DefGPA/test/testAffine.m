
clear all
close all
clc

figDir = '/home/fang/WorkSpace/Clermont-Ferrand/WappredProcrustes/Slides/GPA-Datasets/figures/';


if(0)

dirPath = '../../dataset/Karim_v2/DeformedMesh';


GPA = DGPA_Affine;

GPA.readDir (dirPath);

tic

GPA.run ();

totalUsedTime = toc

%WP.writeShape ('optReferenceAffine');

end


if(1)

load ('../../dataset/face2D/ShapeDataFace.mat')

GPA = DGPA_Affine;

for ii = 1 : size(D,2)
    
    GPA.addPointCloud (D{ii}, diag(V{ii})');
    
end

tic

GPA.run ();

totalUsedTime = toc

GPA.rsdError();

rigidity = GPA.rigidityScore()';
fprintf ( '& %0.1f%% ' , 100*rigidity)
disp('');

GPA


fig = GPA.plotPointCloud ();
fig = tightfig(fig);
%saveas (fig, [figDir, 'affine-face2d-datashape-landmarks'], 'pdf');


fig = GPA.plotRefShape ();
fig = tightfig(fig);
%saveas (fig, [figDir, 'affine-face2d-optimized-shapes'], 'pdf');



end


return




image_origin = zeros(2, GPA.numPointClouds);
max_xy = zeros(2, GPA.numPointClouds);

for ii = 1 : GPA.numPointClouds

    offset = [300; 300];
    f = @(p)  GPA.transformPoints ( p, ii) + offset;
    
    ImgA = imread(['dataset/face2D/', framelist{ii}]);
    
    [B, iorin] = ImageWarp (ImgA, f);
    
    image(ii).B = B;
    image(ii).iorin = iorin;
    
    image_origin(:,ii) = iorin;
    max_xy(:,ii) = [size(B, 2); size(B, 1)];
    
end

% for ii = 1 : GPA.numPointClouds
%     
%     orin = image(ii).iorin;
%     t = image(ii).iorin - image_origin;
%     
%     f = @(p) (t + p);
%     
%     [B, ~] = ImageWarp (image(ii).B, f);
%     
%     figure(ii)
%     imshow(B)
%     
% end

    


xy = min(max_xy, [], 2);

for ii = 1 : GPA.numPointClouds
    
    B = image(ii).B;
    CB{ii} = B(1:xy(2), 1:xy(1), :);
    
end


fig = figure(50);
fig.Position = [50, 800, 1800, 800];
for ii = 1 : GPA.numPointClouds
	subplot(2, ceil(GPA.numPointClouds/2),ii);
    imshow(CB{ii});    
end
fig = tightfig(fig);
%saveas (fig, [figDir, 'affine-face2d-transformed-data'], 'pdf');





eImg = CB{1};
[row, col, cc]= size(eImg);
vv = zeros(1, GPA.numPointClouds);

refImg = CB{1};


for cc = 1 : cc
    for ii = 1 : row
        for jj = 1 : col
            
            for ss = 1 : GPA.numPointClouds
                vv(ss) = CB{ss}(ii,jj,cc);
            end
            eImg(ii, jj, cc) = std(double(vv));
            refImg(ii, jj, cc) = median(vv);
        end
    end
end

fig = figure(100);
fig.Position = [50, 800, 1800, 800];
subplot(2, 3, 1)
imshow(eImg);
title('double')

subplot(2, 3, 2)
imshow(uint8(eImg));
title('uint8')

subplot(2, 3, 3)
imshow(rescale(eImg))
title('rescale')

subplot(2, 3, 4)
imshow(eImg(:,:,1));
title('double - 1st channel')

subplot(2, 3, 5)
imshow(uint8(eImg(:,:,2)));
title('uint8 - 2nd channel')

subplot(2, 3, 6)
imshow(rescale(eImg(:,:,3)));
title('rescale - 3nd channel')

fig = tightfig(fig);
%saveas (fig, [figDir, 'affine-face2d-pixelwise-std'], 'pdf');



fig = figure(200)
imshow(uint8(refImg));
title('Reference image by median');
fig = tightfig(fig);
%saveas (fig, [figDir, 'affine-face2d-reference-shape-Image'], 'pdf');

