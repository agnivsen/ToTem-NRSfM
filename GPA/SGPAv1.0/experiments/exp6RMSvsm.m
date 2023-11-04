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

%Experiment 6 of the paper. Similarity registration: RMS error vs m

addpath('../')
addpath('../util');
options.title='exp6RMSvsm';
options.n=5;
options.d=3;
options.m=50;
options.tau=0.5;
options.sigma=sqrt(0.1);
options.varpar='m';
options.parvalues=round(linspace(9,50,10));
options.iterations=100;
options.methods={'EUC-UPG','EUC-ALL','EUC-ALT'};
genExperiment(options);



