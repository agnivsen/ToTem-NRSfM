% DIFFAROT --- Computes the derivative of a rotation matrix in terms of its
% minimal parameterization for d=2,3.
%
% dif=diffARot(vector,d,similarity)
%
% OUTPUTS
%
%   dif => is the derivative of each term in the rotation matrix, stacked in a
% single column vector, w.r.t. the euler angles of the matrix. 
% Example: d=2 R=[cos(d(1)),-sin(d(1));sin(d(1)),cos(d(1))]; 
%              dif=[-sin(d(1));cos(d(1));-cos(d(1));-sin(d(1))];
%
% INPUTS
%   vector => is the d(d-1)/2 (euclidean) and d(d-1)/2+1 (scaled)
%   parameterization of the rotation matrix.
%   d => is the dimension d=[2,3];
%   similarity={'on','off'} => activate scaled rotation for similarities.
%

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

function dif=diffARot(vector,d,similarity)
if nargin<3
    similarity='off';
end
if(d<2 || d>3)
    disp('Error: Only d=2,3 allowed');
    dif=[];
end

switch(similarity)
    case 'off'
        dvec=length(vector);
        dif=zeros(d*d,dvec);
        switch(d)
            case 2
                dif=[-sin(vector);cos(vector);-cos(vector);-sin(vector)];
            case 3
                [R,dR] = rpyMat (vector);
                dif(:,1)=dR(:,1);
                dif(:,2)=dR(:,2);
                dif(:,3)=dR(:,3);
        end
    case 'on'
        switch(d)
            case 2
                rot=[cos(vector(2:end));sin(vector(2:end));-sin(vector(2:end));cos(vector(2:end))];
                dif=[-sin(vector(2:end));cos(vector(2:end));-cos(vector(2:end));-sin(vector(2:end))];
                dif=[rot,dif];
            case 3
                [R,dR] = rpyMat (vector(2:end));
                dif(:,1)=dR(:,1);
                dif(:,2)=dR(:,2);
                dif(:,3)=dR(:,3);
                dif=[R(:),dif];
        end
end