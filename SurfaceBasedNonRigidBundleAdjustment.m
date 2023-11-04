%% Surface based non-rigid bundle adjustment
% This code replicates our results from figures 8, 9, 10 ,11, and 12 of our ToTem NRSfM 
% paper [1] (pre-print provided in the 'Docs' folder). The method
% implemented here is from section 7 of [1]. 
% -----------------------------------------------------------------------------------
% [1]: Sengupta, Agniva and Adrien Bartoli. "ToTem NRSfM: Object-wise Non-Rigid 
% Structure-from-Motion with a Topological Template" International journal of computer 
% vision, accepted Sep, 2023.

clear all; clear vars; close all; 

addpath('BundleAdjustment\');
addpath('Utils\');
addpath('Parameterisation\');
addpath('GPA\');
addpath('Data\');
addpath('Dependencies\CrustOpen\');

TOPOLOGY_PLANE_ = 0;
TOPOLOGY_CYLINDER_ = 1;
TOPOLOGY_SPHERE_ = 2;

TOPOLOGY = TOPOLOGY_PLANE_;

NUM_TRAINING_FEATURES = 60;
NUM_TRAINING_IMAGES = 11;

dataPath = './Data/';

if(TOPOLOGY == TOPOLOGY_PLANE_)
    dataPath = [dataPath 'SyntheticPlanar/'];
elseif(TOPOLOGY == TOPOLOGY_CYLINDER_)
    dataPath = [dataPath 'SyntheticCylindrical/'];
elseif(TOPOLOGY == TOPOLOGY_SPHERE_)
    dataPath = [dataPath 'SyntheticSpherical/'];
    addpath('./Dependencies/');
    addpath('./Dependencies/gloptipoly3/');
    addpath('./Dependencies/SeDuMi_1_3/');
end

dataWeightInitParam = 0.9; normalWeightInitParam = 0.0008;
reprojectionWeightBA = 0.999; isometryWeightBA =    0.00099;
interpolationDensity = 2000; maxIterationInitParam = 30;
baWeights = [reprojectionWeightBA isometryWeightBA (1 - (reprojectionWeightBA + isometryWeightBA))];
shouldDisplayBAPlots = false;
nK = 7;

synth = SyntheticDataFormatter(dataPath);
synth = synth.setMaxImages(NUM_TRAINING_IMAGES);
synth = synth.setMaxTrainingFeatures(NUM_TRAINING_FEATURES);
synth = synth.setMaxTestFeatures(1000000000);
[DataSubsampled, TestData, K] = synth.simulateSyntheticData();

nFiles = 3;%size(DataSubsampled.p,2);
nPts = size(DataSubsampled.p(1).p,2);

[nng] = getNeighborhoodDuplicated(DataSubsampled,nK);


for iFile = 1:nFiles
    [~, T, cleaned_reconstruction, warpInfo] = InitialRefinement(DataSubsampled.Pgth(iFile).P.', DataSubsampled.v(iFile,:),...
        maxIterationInitParam, TOPOLOGY, dataWeightInitParam, normalWeightInitParam, interpolationDensity);
    
    baData(iFile).template = full(warpInfo.template);
    baData(iFile).target = full(cleaned_reconstruction(:,1:3));
    baData(iFile).T = T;
end

[sparseReconstruction, denseReconstruction, template] = BundleAdjustment(baData, DataSubsampled,...
K, TOPOLOGY, nng, baWeights, shouldDisplayBAPlots, interpolationDensity);

%% Visualising results:

figure(1);
scatter(template(:,1), template(:,2), 'ko', 'filled'); hold on;
for j = 1:nPts
    nN = nng(j,:);
    for iNeighbor = 1:numel(nN)
        q = nN(iNeighbor);
        plot([template(j,1) template(q,1)].',[template(j,2) template(q,2)].','r-');
    end
end
hold on; axis on; grid on; axis equal;
title('Template~$\tau$','Interpreter','latex','FontSize',10);  hold off; pause(0.1);

for iN = 1:nFiles
    sprs = sparseReconstruction(iN).points;
    dns = denseReconstruction(iN).points;
    
    figure(2+iN);
    subplot(1,2,1);
    p1 = plot3(sprs(:,1),sprs(:,2),sprs(:,3), 'o');
    hold on; axis on; grid on; axis equal;
    p1.MarkerFaceColor = [0.9608 0.0941 0.4706];
    p1.MarkerEdgeColor = [0.9608 0.0941 0.4706];
    xlabel('X', 'Interpreter','latex','FontSize',9);
    ylabel('Y', 'Interpreter','latex','FontSize',9);
    zlabel('Z', 'Interpreter','latex','FontSize',9);
    title(['Sparse reconstruction~$\{\mathbf{Q}_{i,j}\}, i = ', num2str(iN),', j \in [1,',num2str(nPts),']$'],'Interpreter','latex','FontSize',9);
    hold off; pause(0.1);
    
    subplot(1,2,2);
    stopDisplayingNewPlot = true;
    pSurf = PlotSurfaceReconstruction(dns);
    pSurf.startShowingAxes();
    pSurf = pSurf.addLightToScene([0 0 -1]);
    pSurf.setMaterial(1);
    pSurf.plotReconstructionNow(stopDisplayingNewPlot);
    pause(0.1); hold off;
    title(['Reconstructed surface~$\varphi_i(\tau), i = ', num2str(iN),'$'],'Interpreter','latex','FontSize',10);
    hold off;
    pause(0.1);
end

