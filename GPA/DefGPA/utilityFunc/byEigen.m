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


% K, the matrix to decompose,
% d, the number of eigen vectors required.
function [V, eigValEigen, eigenValueAll] = byEigen (K, d, str)

offset = 0;

[Q, D ] = eig( (K+K')/2 );

eigVal = diag(D);

[eigVal, Index] = sort(eigVal, 'descend');  Q = Q (:, Index);

if  strcmp (str, 'largest')
    
    V = Q (:, 1+offset:d+offset);
    
    eigValEigen =  eigVal (1+offset:d+offset);
    
end

if strcmp (str, 'smallest')
    
    V = Q (:, end-d+1-offset:end-offset);
    
    eigValEigen =  eigVal (end-d+1-offset:end-offset);
    
end

eigenValueAll = eigVal;

%             if  norm( V' * V - eye(d), 'fro' ) > 1e-12
%
%                 fprintf(2, 'Eigen vectors are not orthnormal!.');
%
%             end

end