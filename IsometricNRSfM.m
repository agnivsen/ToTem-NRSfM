%% Unconstrained NRSfM
% This code runs unconstrained NRSfM as described in Sec. 4 of our ToTem
% NRSfM paper [1] (pre-print provided in the 'Docs' folder). This script
% exactly replicates our results in the first two columns of Table 3 in
% [1]. Our method, requires as input, not only the point correspondences
% and intrinsics, but also an initialising depth, since our proposed method
% is a non-convex optimisation. As detailed in Sec. 8.1 of [1], we obtain
% these initialising depth from two sources, [2] and [3]. These initial
% depths are provided along with the data in the 'Data' folder. To run this
% code in your own data, please obtain the source code of [2] or [3] and
% use them to initialise our method. Alternatively, you can skip providing
% any initial depth to our method and the code shall initialise all depths
% to unity; this shall produce a feasible solution but the accuracy may
% suffer.
%
% This code has been tested on Ubuntu 20.04 (Focal Fossa) with Matlab 2022b
% and Matlab 2016b and on Windows 11 with Matlab 2022a. 
% -----------------------------------------------------------------------------------
% [1]: Sengupta, Agniva and Adrien Bartoli. "ToTem NRSfM: Object-wise Non-Rigid 
% Structure-from-Motion with a Topological Template" International journal of computer 
% vision, accepted Sep, 2023.
%
% [2]: Chhatkuli, A., Pizarro, D., Collins, T., & Bartoli, A. (2017). Inextensible non-rigid 
% structure-from-motion by second-order cone programming. IEEE transactions on pattern 
% analysis and machine intelligence, 40(10), 2428-2441.
%
% [3]: Ji, P., Li, H., Dai, Y., & Reid, I. (2017). " Maximizing rigidity" revisited: a convex 
% programming approach for generic 3D shape reconstruction from multiple perspective 
% views. In Proceedings of the IEEE International Conference on Computer Vision 


clear all; close all; clear vars;

% adding data path
addpath('Data/');
addpath('NRSfM/');
addpath('Utils/');

% Four options for data:
% data options: 'hulk', 'whiteCartoonTshirt', 'kinectPaper', 'cushion'
dataName = 'cushion';

% loading data and intrinsics
data = load(['Data/data_' dataName '.mat']);
K = load(['Data/intrinsics_' dataName '.txt']);

% loading initial depths from 
initFiles = load(['Data/' dataName '_Initialisation.mat']);

% Invoking NRSfM
nrMDH = NRSfM(data); % setting the data
nrMDH.setIntrinsics(K); % setting the intrinsics
 % initial depth from [2]
nrMDH.setInitialDepth(initFiles.initialisations.MDH_SOCP.initialDepth);
% invoking NRSfM computation, argument is the no. of neighbours
nrMDH.computeNRSfM(initFiles.initialisations.MDH_SOCP.nK); 
% fetching the reconstructed 3D points
[reconstructionOurs_MDH] = nrMDH.get3DReconstruction(); 

% Invoking NRSfM, yet again
nrMLH = NRSfM(data); % setting the data
nrMLH.setIntrinsics(K); % setting the intrinsics
% initial depth from [3]
nrMLH.setInitialDepth(initFiles.initialisations.MLH_SDP.initialDepth); 
% invoking NRSfM computation, argument is the no. of neighbours
nrMLH.computeNRSfM(initFiles.initialisations.MLH_SDP.nK);
% fetching the reconstructed 3D points
[reconstructionOurs_MLH] = nrMLH.get3DReconstruction();



%% Plotting results:
% ===> Requires Matlab 2019b or later
figure;
tPlot = tiledlayout('flow');
for n = 1:size(data.Pgth,2)
    R1 = reconstructionOurs_MDH{n};
    R2 = reconstructionOurs_MLH{n};
    GT = data.Pgth(n).P.';
    nexttile;
    plot3(R1(:,1), R1(:,2), R1(:,3), 'ro', 'MarkerFaceColor','r'); hold on;
    plot3(R2(:,1), R2(:,2), R2(:,3), 'bo', 'MarkerFaceColor','b');
    plot3(GT(:,1), GT(:,2), GT(:,3), 'ko', 'MarkerFaceColor','k');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis on; grid on; axis equal;
    title(['$n =~$' num2str(n)], 'Interpreter','latex', 'FontSize',7);
    hold off; pause(0.1);
end

leg = legend({'uNRS\textit{f}M-MDH-SOCP', 'uNRS\textit{f}M-MLH-SDP', 'GT'}, ...
    "Interpreter","latex", "FontSize",9);
leg.Layout.Tile = 'north';
leg.Orientation = 'horizontal';

title(tPlot, ['Data = ' dataName], 'Interpreter', 'Latex');
