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


% a function used to construct a sparse matrix by given sublocks
% this function put a matrix A at the position (row, col) of a sparse matrix.
% use " sparse (i, j, v) " to check the resulting sparse matrix
function [i, j, v] = putMatrixBlock (A, row, col)

[m, n] = size(A);

v = full (reshape (A, 1, m * n));

i = kron (ones(1,n), row:row+m-1);

j = kron (col:col+n-1, ones(1,m));

end
