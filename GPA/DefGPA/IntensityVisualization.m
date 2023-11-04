% Copyright 2021 Fang Bai <fang dot bai at yahoo dot com>
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.


function [tformImgCell, varImg, refImg] = IntensityVisualization (DGPA_Handle, ImgFileCell, GlobalRotation, OriginOffset)

close all;

image_origin = zeros(2, DGPA_Handle.numPointClouds);

max_xy = zeros(2, DGPA_Handle.numPointClouds);

for ii = 1 : DGPA_Handle.numPointClouds
    
    f = @(p)  GlobalRotation * DGPA_Handle.transformPoints ( p, ii) + [OriginOffset(1); OriginOffset(2)];
    
    ImgA = imread(ImgFileCell{ii});
    
    [B, iorin] = ImageWarp (ImgA, f);
    
    image(ii).B = B;
    image(ii).iorin = iorin;
    
    image_origin(:,ii) = iorin;
    max_xy(:,ii) = [size(B, 2); size(B, 1)];
    
end


xy = min(max_xy, [], 2);

for ii = 1 : DGPA_Handle.numPointClouds
    B = image(ii).B;
    tformImgCell{ii} = B(1:xy(2), 1:xy(1), :);
end


varImg = tformImgCell{1};
[row, col, cc]= size(varImg);
vv = zeros(1, DGPA_Handle.numPointClouds);
refImg = tformImgCell{1};
for cc = 1 : cc
    for ii = 1 : row
        for jj = 1 : col
            for ss = 1 : DGPA_Handle.numPointClouds
                vv(ss) = tformImgCell{ss}(ii,jj,cc);
            end
            varImg(ii, jj, cc) = std(double(vv));
            refImg(ii, jj, cc) = median(vv);
        end
    end
end

end