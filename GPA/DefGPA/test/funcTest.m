clear allc
close all
clc

img = imread ( 'dataset/face2D/deform0090.jpg' );



m = [ 2, 3, 100;
        1, 4, 200;
        0, 0, 1];
    

tform = affine2d(m');
C = imwarp(img, tform);


figure(3)
imshow(C)





[B, origin] = ImageWarp (img, @warpFunc);

figure(4)
imshow(B)

origin

c1 =  [2 3; 1 4] * [1; 1]
c2 =  [2 3; 1 4] * [720; 1]
c3 =  [2 3; 1 4] * [1; 576]
c4 =  [2 3; 1 4] * [720; 576]



function p =  warpFunc (pp)
    
    %p = [2 -1; 1 4] * pp;
    
    p = [2 3; 1 4] * pp;
    
end