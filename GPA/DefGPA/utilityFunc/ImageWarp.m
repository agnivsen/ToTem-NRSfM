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


function [B, origin] = ImageWarp (A, warpFunc)

[row, col, cc] = size(A);

[X, Y ] = meshgrid(1:col, 1:row);

for ii = 1 : row
    
    for jj = 1 : col
        
        po = warpFunc ([X(ii,jj); Y(ii,jj)]);
        
        X(ii,jj) = po(1);
        
        Y(ii,jj) = po(2);
        
    end
    
end

max_x = max(max(X));
max_y = max(max(Y));

min_x = min(min(X));
min_y = min(min(Y));

% set to origin if possible
min_x = 1;
min_y = 1;

[Xr, Yr] = meshgrid(min_x:max_x, min_y:max_y);

for ii = 1 : cc
    
    V = griddata(X, Y, double(A(:,:,ii)), Xr, Yr);
    
    B(:,:,ii) = uint8(V);
    
end

origin = [min_x; min_y];

end
