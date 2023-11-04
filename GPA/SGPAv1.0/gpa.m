% gpa --- Generalized Procrustes Analysis with affine or similarity/euclidean transformations.
%
% This function outputs the solution of the GPA. It implements classical
% alternation method and the stratified solution of the generalized
% procrustes anaylisis proposed by Bartoli,Pizarro,Loog paper "Stratified
% Generalized Procrustes Analysis" BMVC2010.
%
% T=gpa(D,V,options)
%
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
% D is the shape cell array.
%       D{i} is the dxm shape matrix of shape i. (D{i}=[Di1,....,Dim], Dij=[v1,..,vd]')
% V is the matrix with missing data V(i,j)=1 if point j at shape i is
% missing.
% options is an input structure with the following fields.
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

function T=gpa(D,V,options)

if(nargin==1)
    [D,V,options,n,d,m]=ProcessArgs(D);
elseif(nargin==2)
    [D,V,options,n,d,m]=ProcessArgs(D,V);
else
    [D,V,options,n,d,m]=ProcessArgs(D,V,options);
end

if(n<2)
    disp('Error: n>1 shapes are needed...');
    T=[];
    return
end

switch(lower(options.method))
    case 'aff-fct'
        opti.Ref='ni';
        opti.Init='fct';
        opti.verbose=options.verbose;
        opti.maxiter=options.maxiter;
        T=AffineRegistration(D,V,opti);
    case 'aff-ref'
        opti.Ref='ni';
        opti.Init='ref';
        opti.verbose=options.verbose;
        opti.maxiter=options.maxiter;
        T=AffineRegistration(D,V,opti);
    case 'aff-all'
        opti.Ref='lm';
        opti.Init='ref';
        opti.verbose=options.verbose;
        opti.maxiter=options.maxiter;
        T=AffineRegistration(D,V,opti);
    case 'sim-upg'
        opti.Ref='ni';
        opti.upgrade='fromreg';
        opti.similarity='on';
        opti.method='stratified';
        opti.verbose=options.verbose;
        opti.maxiter=options.maxiter;
        T=EuclideanRegistration(D,V,opti);
    case 'sim-all'
        opti.Ref='lm';
        opti.upgrade='fromreg';
        opti.similarity='on';
        opti.method='stratified';
        opti.verbose=options.verbose;
        opti.maxiter=options.maxiter;
        T=EuclideanRegistration(D,V,opti);
    case 'sim-alt'
        opti.similarity='on';
        opti.method='alternation';
        opti.verbose=options.verbose;
        opti.maxiter=options.maxiter;
        T=EuclideanRegistration(D,V,opti);
    case 'euc-upg'
        opti.Ref='ni';
        opti.upgrade='fromreg';
        opti.similarity='off';
        opti.method='stratified';
        opti.verbose=options.verbose;
        opti.maxiter=options.maxiter;
        T=EuclideanRegistration(D,V,opti);
    case 'euc-all'
        opti.Ref='lm';
        opti.upgrade='fromreg';
        opti.similarity='off';
        opti.method='stratified';
        opti.verbose=options.verbose;
        opti.maxiter=options.maxiter;
        T=EuclideanRegistration(D,V,opti);
    case 'euc-alt'
        opti.similarity='off';
        opti.method='alternation';
        opti.verbose=options.verbose;
        opti.maxiter=options.maxiter;
        T=EuclideanRegistration(D,V,opti);
end

end

function [D,V,options,n,d,m]=ProcessArgs(D,V,options)

n=length(D);
if(n>0)
    [d,m]=size(D{1});
end
if(nargin<3)
    options.method='SIM-ALL';
    options.maxiter=200;
    options.verbose='off';
    if(nargin<2)
        for i=1:n
            V{i}=ones(1,m);
        end
    end
else

    if(~isfield(options,'method'))
        options.method='SIM-ALL';
    end
    if(~isfield(options,'maxiter'))
        options.maxiter=200;
    end
    if(~isfield(options,'verbose'))
        options.verbose='off';
    end
end


for i = 1:n
    [d1, d2] = size(V{i});
    if(d1 == d2)
        V{i} = diag(V{i})';
    end
end

end

