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


% Run experiments with synthetic data. (The experiments are very time consuming !.)
% For a quick demo try to reduce the number of iterations on each file
%Experiment 1          (Affine registration with rigid shapes)
exp1RMSvsd
%%%UNCOMMENT TO RUN ALL EXPERIMENTS
% exp1RMSvsm
% exp1RMSvsn
% exp1RMSvssigma
% % Experiment 2       (Affine registration with non-rigid shapes)
% exp2RMSvsd
% exp2RMSvsm
% exp2RMSvsn
% exp2RMSvssigma
% % Experiment 3
% exp3nonrigidRMSvstau (Affine registration with non-rigid shapes and missing data)
% exp3rigidRMSvstau    (Affine registration with rigid shapes and missing data)
% % Experiment 4       (Similarity registration with rigid shapes and missing data)
% exp4RMSvsd
% exp4RMSvsm
% exp4RMSvsn
% exp4RMSvssigma
% exp4RMSvstau
% %Experiment 5        (Similarity registration with rigid shapes and missing data)
% exp5RMSvsd
% exp5RMSvsm
% exp5RMSvsn
% exp5RMSvssigma
% exp5RMSvstau
% %Experiment 6        (Similarity registration with non-rigid shapes and missing data)
% exp6RMSvsd
% exp6RMSvsm
% exp6RMSvsn
% exp6RMSvssigma
% exp6RMSvstau
% %Experiment 7        (Euclidean registration with non-rigid shapes and missing data)
% exp6RMSvsd
% exp6RMSvsm
% exp6RMSvsn
% exp6RMSvssigma
% exp6RMSvstau

%% Plot a specific experiment (eg exp1RMSvsd) in figure 1 and saved as
%% exp1RMSvsd.fig

nameexp1RMSvsd=''; % This will show an open file dialog to select mat file in mats/ directory
                   % If a filename is given the function will not show the
                   % dialog (eg
                   % nameexp1RMSvsd='exp1RMSvsd-2010-11-23-8-46-46.mat'
plotExperiment(nameexp1RMSvsd,1,'exp1RMSvsd');


