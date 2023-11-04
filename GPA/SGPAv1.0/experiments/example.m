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

clear all;
addpath('../util');
addpath('../');
% Usage Example
% n=number of shapes
% m=number of d-dimensional points inside each shape
% d=dimension of the points of each shape (e.g. d=2, the points are two-dimensional)
%Each shape is inside a cell array D where D{i} is a dxm matrix
%The missing data is included in the array nxm array V
%where V(i,j)=1 if the point is visible and V(i,j)=0 if is hidden.

n=50;
m=40;
d=3;
for i=[1:n]
D{i}=randn(d,m);
end
V=ones(n,m);

%options is an input structure with the following fields.
%      options.method={
%           'AFF-FCT' => Affine registration with closed-form factorization method (no
%                      missing data allowed)
%           'AFF-REF' => Affine registration with closed-form in
%                       reference-space. (missing data is allowed).
%           'AFF-ALL' => Affine registration with iterative refinement.
%                      (AFF-REF is used as initialization except when no missing data
%                      is found).
%           'SIM-UPG' => Similarity registration using upgrade from registration parameters 
%                      of AFF-ALL.
%           'SIM-ALL' => SIM-UPG + iterative refinement in data-space.
%           (DEFAULT)
%           'SIM-ALT' => Classical Alternation method for similarity transformations.
%           'EUC-UPG' => Euclidean registration using upgrade from registration parameters 
%                      of AFF-ALL.
%           'EUC-ALL' => EUC-UPG + iterative refinement in data-space.
%           'EUC-ALT' => Classical Alternation method for Euclidean case
%           transformations.
%           }
%      options.verbose={'on','off'}
%       'on' => show debug information and error per iteration
%       'off'=> don't show debug information.
%      options.maxiter => maximum number of iterations (=200
%      DEFAULT)

options.method='AFF-ALL';
options.verbose='on';
T=gpa(D,V,options);

% T is a structure with fields:
%       T.A is a cell-array with the rotational part of the transformations
%       T.a is a cell-array with the translational part
%       T.S is the reference shape.
%       T.errorvec is the data-space error obtained across iterations
%       T.niter is the number of iterations required.
%       T.time is the processing time of the algorithm.
%           T.time(1) is the total ellapsed time.
%           T.time(2) is the ellapsed time in the Initialization.
%           T.time(3) is the average time per iteration.
%           T.time(4) is the total time ellapsed in the refinement.

