clear all
close all
clc


addpath('../DefGPA', '../');
addpath('../SGPAv1.0', '../SGPAv1.0/util');

fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end


figformat = 'pdf';


datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug'}

centerePerAxis = [3, 5, 7];

dimFeatureSpace = {centerePerAxis.^2,  centerePerAxis.^2, centerePerAxis.^2, centerePerAxis.^3, centerePerAxis.^3, 2*centerePerAxis.^2};


sp1 = [0.5, 1, 2, 5, 10, 25, 50, 75, 100, 500, 1000, 10000, 100000];
sp2 = [1e-5, 1e-4, 1e-3, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.1, 1];
sp3 = [1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10, 50, 100 ];

sp4 = [1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10, 50, 100 ];

sp5 = [1e-3, 2.5e-3, 5e-3, 7.5e-3, 1e-2, 2.5e-2, 5e-2, 7.5e-2, 0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10, 50, 100 ];

sp5 = [1e-2, 2.5e-2, 5e-2, 7.5e-2, 0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10, 50, 100 ];

smoothParamCell = {sp1, sp1, sp1, sp5, sp3, sp4};

bestParams = {10, 100, 100, 0.1, 0.01, 0.01};

ticks1 = [1, 1e1, 1e2, 1e3, 1e4, 1e5];
ticks2 = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
ticks3 = [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2];

ticks4 = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2];

ticksCell ={ticks1, ticks1, ticks1, ticks3, ticks3, ticks4};

n = [1,2,3,4, 5, 6]

n = [1, 2, 3, 4];

n = 1:2;




ColorMatrix =  [0    0.4470    0.7410;
                            0.8500    0.3250    0.0980;
                            0.9290    0.6940    0.1250;
                            0.4940    0.1840    0.5560;
                            0.4660    0.6740    0.1880;
                            0.3010    0.7450    0.9330;
                            0.6350    0.0780    0.1840;];

newColorMatrix = [0.83 0.14 0.14;
                        1.00 0.54 0.00;
                        0.47 0.25 0.80;
                        0.25 0.80 0.54];

newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};                    


ColorMatrix = jet(4);


usedColor1 = ColorMatrix(1,:);
usedColor2 = ColorMatrix(2,:);
usedColor3 = ColorMatrix(3,:);
usedColor4 = ColorMatrix(4,:);
% usedColor4 = 'g';


% usedColor1 = newcolors{1};
% usedColor2 = newcolors{5};
% usedColor3 = newcolors{7};
% usedColor4 = 'g';



for dtc = 1:6
    
datasetName = datasetCell{dtc}

smoothParam = smoothParamCell{dtc};

xticksval = ticksCell{dtc};


[DataCell, VisVecCell] = readDataset (datasetName);


if dtc == 5
    if (1)
        [RMSE_r, RMSE_d, RigScore, crossValidationError] = runTPS_Params (DataCell, VisVecCell, dimFeatureSpace{dtc}, smoothParam);
        save('Liver_TPS_Params_RMSE_CVE.mat', 'RMSE_r', 'RMSE_d', 'RigScore', 'crossValidationError');
    else
        load('Liver_TPS_Params_RMSE_CVE.mat', 'RMSE_r', 'RMSE_d', 'RigScore', 'crossValidationError');
    end
else
    [RMSE_r, RMSE_d, RigScore, crossValidationError] = runTPS_Params (DataCell, VisVecCell, dimFeatureSpace{dtc}, smoothParam);
end
% if dtc == 6
%     %save('ToyRig_TPS_Params_RMSE_CVE.mat', 'RMSE_r', 'RMSE_d', 'RigScore', 'crossValidationError');
%     load('ToyRig_TPS_Params_RMSE_CVE.mat', 'RMSE_r', 'RMSE_d', 'RigScore', 'crossValidationError');
% end


% AFFINE
GPA_AFFINE = DefGPA('AFFINE');
GPA_AFFINE.addDataByArray(DataCell, VisVecCell);
GPA_AFFINE.run ();
GPA_AFFINE.rsdError();
GPA_AFFINE.rigidityScore();
AFF_r.RMSE_r = GPA_AFFINE.rsd_ref;
AFF_r.RMSE_d = GPA_AFFINE.rsd_dat;
AFF_r.rigidity = mean(GPA_AFFINE.rigidity);
AFF_r.crossValidationError  = computeCrossValidationError(DataCell, VisVecCell, 0, inf, GPA_AFFINE.mShape);


close all
hFig = figure(dtc);
set(hFig, 'Position', [100, 50, 500, 200]);
%colorOptions = {'r', 'b', 'm', 'g'};
colorOptions = {usedColor1, usedColor2, usedColor3, usedColor4};
lineStyleOption = {':', '-', ':', '-'};
markerOption = {'s', '.', '.'};
lineWidth = 1.5;


xval = smoothParam;
subplot(1, 2, 1);
hold on
ahd = plot (xval, AFF_r.RMSE_d * ones(1, max(size(xval))), 'LineStyle', lineStyleOption{1}, 'LineWidth', lineWidth, 'Color', colorOptions{max(size(centerePerAxis))+1}, 'Marker', markerOption{1}, 'MarkerSize', 8);
ahr = plot (xval, AFF_r.RMSE_r * ones(1, max(size(xval))), 'LineStyle', lineStyleOption{2}, 'LineWidth', lineWidth, 'Color', colorOptions{max(size(centerePerAxis))+1}, 'Marker', markerOption{2}, 'MarkerSize', 12);
for ii = 1 : max(size(centerePerAxis))
    hd{ii} = plot(xval, RMSE_d(ii,:), 'LineStyle', lineStyleOption{1}, 'LineWidth', lineWidth, 'Color', colorOptions{ii}, 'Marker', markerOption{1}, 'MarkerSize', 8);
end
for ii = 1 : max(size(centerePerAxis))
    hr{ii} = plot(xval, RMSE_r(ii,:), 'LineStyle', lineStyleOption{2}, 'LineWidth', lineWidth, 'Color', colorOptions{ii}, 'Marker', markerOption{2}, 'MarkerSize', 12);
end
ahx = xline (bestParams{dtc}, '-', {['\theta = ', num2str(bestParams{dtc})]}, 'LineWidth', lineWidth, 'FontSize', 14, 'LabelVerticalAlignment', 'bottom');
hold off

set(gca, 'XScale', 'log');
xticks(xticksval);
grid on;

%xlabel('$\theta$', 'FontSize', 10, 'Interpreter', 'latex');
ylabel('RMSE (pixels)', 'FontSize', 10);
set(gca, 'FontSize', 10);
box on



subplot(1, 2, 2);
hold on
ah = plot (xval, AFF_r.crossValidationError * ones(1, max(size(xval))), 'LineStyle', lineStyleOption{2}, 'LineWidth', lineWidth, 'Color', colorOptions{max(size(centerePerAxis))+1}, 'Marker', markerOption{2}, 'MarkerSize', 12);
for ii = 1 : max(size(centerePerAxis))
    h{ii} = plot(xval, crossValidationError(ii,:), 'LineStyle', lineStyleOption{2}, 'LineWidth', lineWidth, 'Color', colorOptions{ii}, 'Marker', markerOption{2}, 'MarkerSize', 12);
end
ahx = xline (bestParams{dtc}, '-', {['\theta = ', num2str(bestParams{dtc})]}, 'LineWidth', lineWidth, 'FontSize', 14, 'LabelVerticalAlignment', 'bottom');
hold off

set(gca, 'XScale', 'log');
xticks(xticksval);
grid on;

%xlabel('$\theta$', 'FontSize', 10, 'Interpreter', 'latex');
ylabel('CVE (pixels)', 'FontSize', 10);
set(gca, 'FontSize', 10);
box on



hFig = tightfig(hFig);
fig_saved_name = [fdir, 'exp_TPS_Smooth_Parms_' datasetName];
%saveas (hFig, fig_saved_name, figformat );
print(hFig,  '-dpdf', fig_saved_name, '-painters');


end



if(0)
    
fhandles = [ahd, ahr];
legendNames = {'AFF-DAT', 'AFF-REF'};
for ii = 1 : max(size(centerePerAxis))
    fhandles(end+1) = hd{ii};
    legendNames{end+1} = ['TPS(', num2str(centerePerAxis(ii)), ')-DAT'];
    fhandles(end+1) = hr{ii};
    legendNames{end+1} = ['TPS(', num2str(centerePerAxis(ii)), ')-REF'];
end
hL1 = legend(fhandles, legendNames, 'FontSize', 8, 'Orientation', 'horizontal', 'Box', 'on', 'LineWidth', 1);



llength = 0.80;
lwidth = 0.1;
legPos = [0.51-llength/2, -0.10-lwidth, llength, lwidth];
set(hL1,'Position', legPos,'Units', 'normalized');


set(hFig, 'Position', [100, 50, 960, 200 * (1 + lwidth +0.2) ]);
hFig = tightfig(hFig);
fig_saved_name = [fdir, 'exp_TPS_Smooth_Parms_caption'];
print(hFig,  '-dpdf', fig_saved_name, '-painters');

%%
%pdfcrop exp_TPS_Smooth_Parms_caption.pdf --margins '-5 -205 -5 0'


end








function [RMSE_r, RMSE_d, RigScore, crossValidationError] = runTPS_Params (DataCell, VisVecCell, dimFeatureSpace, smoothParam)

allM = zeros(max(size(dimFeatureSpace)), max(size(smoothParam)));
RMSE_r = allM;
RMSE_d = allM;
RigScore = allM;
crossValidationError = allM;

for ii = 1 : max(size(dimFeatureSpace))
    
    for jj = 1 : max(size(smoothParam))
        
        % TPS warp
        GPA_WARP = DefGPA('TPS');
        GPA_WARP.addDataByArray(DataCell, VisVecCell);
        GPA_WARP.dimFeatureSpace = dimFeatureSpace(ii);
        GPA_WARP.flagJointlyOptimizeLamdaRigidTransformation = false;
        GPA_WARP.smoothParam = smoothParam(jj);
        GPA_WARP.run(smoothParam(jj) );
        GPA_WARP.rsdError();
        GPA_WARP.rigidityScore();
        
        RMSE_r(ii,jj) = GPA_WARP.rsd_ref;
        RMSE_d(ii,jj) = GPA_WARP.rsd_dat;
        RigScore(ii,jj) = mean(GPA_WARP.rigidity);

        crossValidationError(ii, jj)  = computeCrossValidationError(DataCell, VisVecCell, dimFeatureSpace(ii), smoothParam(jj), GPA_WARP.mShape);
                
        clear GPA_WARP
        
    end
    
end

end




