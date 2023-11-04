% DATASPACEERROR --- Computes the residual error of the procrustes solution
% in the data space.
%
% error=dataspaceerror(T,D,V,verbose)
%
% error is the data space scalar error error=1/(n*m) sum_i sum_j || Dij- T_i(Sj)||^2_2
% T is a structure with fields:
%       T.A is a cell-array with the rotational part of the transformations
%       T.a is a cell-array with the translational part
%       T.S is the reference shape.
%   It is compatible with the output structure from gpa function.
% D is the shape cell array.
%       D{i} is the dxm shape matrix of shape i. (D{i}=[Di1,....,Dim], Dij=[v1,..,vd]')
% V is the matrix with missing data V{i}(j,j)=1 if point j at shape i is
% missing.
% verbose=['on','off} to activate warning if some of the transformations
%                     has negative determinant.

% Copyright (C)2010  Adrien Bartoli (1), Daniel Pizarro (2) and Marco Loog
% (3)
% (1) Clermont Université - Clermont-Ferrand, France.
% (2) Universidad de Alcala- Alcala de Henares, Spain
% (3) Delft University - Nederlands
% Send bug reports and feedback to: dani.pizarro@gmail.com
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


function error=dataspaceerror(T,D,V,verbose)
if(nargin<4)
    verbose='off';
end
n=length(D);
[d,m]=size(D{1});
error=0;
alpha=0;
w=sign(det(T.A{1}));
for i=1:n
    if(w~=sign(det(T.A{i})))
        dispv(sprintf('Shape %d has negative determinant',i),verbose);
    end
    for j=1:m
        if(V{i}(j))
            error=error+norm(D{i}(:,j)-T.A{i}*T.S(:,j)-T.a{i})^2;
            alpha=alpha+1;
        end
    end
end
error=error/alpha;