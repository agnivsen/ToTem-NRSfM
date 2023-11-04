% REFRENCESPACEERROR --- Computes the residual error of the procrustes solution in the reference 
%                        space.
%
% error=dataspaceerror(T,D,V)
%
% error is the reference space scalar error error=1/(n*m) sum_i sum_j || Sj-T^{-1}_i(Dij)||^2_2
% T is a structure with fields:
%       T.A is a cell-array with the rotational part of the transformations
%       T.a is a cell-array with the translational part
%       T.S is the reference shape.
%   It is compatible with the output structure from gpa function.
% D is the shape cell array.
%       D{i} is the dxm shape matrix of shape i. (D{i}=[Di1,....,Dim], Dij=[v1,..,vd]')
% V is the matrix with missing data V{i}(j,j)=1 if point j at shape i is
% missing.

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

function [error,T2]=referencespaceerror(T,D,V)
n=length(D);
[d,m]=size(D{1});
error=0;
alpha=0;
M=T.A{1};
T.S=M*T.S;
toff=T.S*ones(m,1)./m;
T.S=T.S-toff*ones(1,m);
T2=T;
for i=[1:n]
    T.A{i}=M*inv(T.A{i});
    T.a{i}=-T.A{i}*T.a{i}-toff;
    T2.A{i}=inv(T.A{i});
    T2.a{i}=-inv(T.A{i})*T.a{i};
end
T2.S=T.S;
for i=1:n
    for j=1:m
        if(V{i}(j))
            error=error+norm(T.S(:,j) - T.A{i}*D{i}(:,j)-T.a{i})^2;
            alpha=alpha+1;
        end
    end

end
error=error/alpha;