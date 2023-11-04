clear all
close all
clc


addpath('../DefGPA', '../');


fdir = '../../PaperDraft/figure/';
if ~exist(fdir, 'dir')
    mkdir(fdir);
end


% ----- Dataset Information ---

datasetCell = {'Face', 'HandBag', 'PillowCover', 'LiTS', 'Liver', 'StereoToyRug', 'TOPACS'};


DefGPA_TPS_bestParams = {10, 100, 100, 0.1, 0.01, 0.01};
dimensions = {2, 2, 2, 3, 3, 3, 3}; 


KGPA_bestQuantilePercentParams = [0.2,  0.2,  0.2,  0.2,  0.2,  0.2];
KGPA_bestSmoothParams = [0.05, 0.05, 0.05, 0.05,  0.05, 0.05];


dtc = 2;
datasetName = datasetCell{dtc};

num_correspondences_to_remove = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600];


TPS_smooth_param = DefGPA_TPS_bestParams{dtc};
TPS_control_points_dim = dimensions{dtc};

Kernel_p_quantile = KGPA_bestQuantilePercentParams(dtc);
Kernel_smooth_param = KGPA_bestSmoothParams(dtc);


IJCV2022_CVE_Leave1 = cell(1, length(num_correspondences_to_remove));
IJCV2022_CVE_Nfold = cell(1, length(num_correspondences_to_remove));

RSS2022_CVE_Leave1 = cell(1, length(num_correspondences_to_remove));
RSS2022_CVE_Nfold = cell(1, length(num_correspondences_to_remove));


% By setting Nfold to infinity, the method
% run_N_Fold_CrossValidation(Nfold)
% generates Leave-1-out cross-validation.
Nfold = 20;

MC_Trials = 20;


for c_cnt = 1 : length(num_correspondences_to_remove)
    
    [DataCellFull, VisVecCellFull] = readDataset (datasetName);
    
    IJCV2022_CVE_Leave1{c_cnt} = [];
    IJCV2022_CVE_Nfold{c_cnt} = [];

    RSS2022_CVE_Leave1{c_cnt} = [];
    RSS2022_CVE_Nfold{c_cnt} = [];
    
    fprintf(1, "process case: %d \n",  num_correspondences_to_remove(c_cnt));

    
    for mc_cnt = 1 : MC_Trials
        
        [DataCell, VisVecCell] = RemoveCorrespondences (DataCellFull, VisVecCellFull, num_correspondences_to_remove(c_cnt));        
        
        % this use the Lambda estimation method in IJCV2022
        flagJointlyOptimizeLamdaRigidTransformation = false;
        
        [cve_rigid_leave1, cve_affine_leave1, cve_tps3_leave1, cve_tps5_leave1, cve_tps7_leave1, cve_kernel_leave1] = runAllMethods (DataCell, VisVecCell, ...
            TPS_smooth_param, TPS_control_points_dim, ...
            Kernel_p_quantile, Kernel_smooth_param, ...
            flagJointlyOptimizeLamdaRigidTransformation, inf);
        
        [cve_rigid_Nfold, cve_affine_Nfold, cve_tps3_Nfold, cve_tps5_Nfold, cve_tps7_Nfold, cve_kernel_Nfold] = runAllMethods (DataCell, VisVecCell, ...
            TPS_smooth_param, TPS_control_points_dim, ...
            Kernel_p_quantile, Kernel_smooth_param, ...
            flagJointlyOptimizeLamdaRigidTransformation, Nfold);
        
        IJCV2022_CVE_Leave1{c_cnt} = [IJCV2022_CVE_Leave1{c_cnt};
            [cve_rigid_leave1, cve_affine_leave1, cve_tps3_leave1, cve_tps5_leave1, cve_tps7_leave1, cve_kernel_leave1] ];
        IJCV2022_CVE_Nfold{c_cnt} = [IJCV2022_CVE_Nfold{c_cnt};
            [cve_rigid_Nfold, cve_affine_Nfold, cve_tps3_Nfold, cve_tps5_Nfold, cve_tps7_Nfold, cve_kernel_Nfold] ];
        
        
        
        % this use the novel Lambda estimation method in RSS2022
        flagJointlyOptimizeLamdaRigidTransformation = true;
        
        [cve_rigid_leave1, cve_affine_leave1, cve_tps3_leave1, cve_tps5_leave1, cve_tps7_leave1, cve_kernel_leave1] = runAllMethods (DataCell, VisVecCell, ...
            TPS_smooth_param, TPS_control_points_dim, ...
            Kernel_p_quantile, Kernel_smooth_param, ...
            flagJointlyOptimizeLamdaRigidTransformation, inf);
        
        [cve_rigid_Nfold, cve_affine_Nfold, cve_tps3_Nfold, cve_tps5_Nfold, cve_tps7_Nfold, cve_kernel_Nfold] = runAllMethods (DataCell, VisVecCell, ...
            TPS_smooth_param, TPS_control_points_dim, ...
            Kernel_p_quantile, Kernel_smooth_param, ...
            flagJointlyOptimizeLamdaRigidTransformation, Nfold);
        
        RSS2022_CVE_Leave1{c_cnt} = [RSS2022_CVE_Leave1{c_cnt};
            [cve_rigid_leave1, cve_affine_leave1, cve_tps3_leave1, cve_tps5_leave1, cve_tps7_leave1, cve_kernel_leave1] ];
        RSS2022_CVE_Nfold{c_cnt} = [RSS2022_CVE_Nfold{c_cnt};
            [cve_rigid_Nfold, cve_affine_Nfold, cve_tps3_Nfold, cve_tps5_Nfold, cve_tps7_Nfold, cve_kernel_Nfold] ];
        
    end

end


save(['Lambda_Methods_by_Missing_Correspondences_', datasetName, '.mat'], "RSS2022_CVE_Leave1", "RSS2022_CVE_Nfold", "IJCV2022_CVE_Leave1", "IJCV2022_CVE_Nfold", "num_correspondences_to_remove");






use_percentage = true;


CVE_Leave1_EUC = [];
CVE_Leave1_AFF = [];
CVE_Leave1_TPS3 = [];
CVE_Leave1_TPS5 = [];
CVE_Leave1_TPS7 = [];
CVE_Leave1_Kernel = [];

for ii = 1 : length(num_correspondences_to_remove)
    ijcv_case = IJCV2022_CVE_Leave1{ii};
    rss_case = RSS2022_CVE_Leave1{ii};    
    err_per_case = (ijcv_case - rss_case) ./ ijcv_case;
    CVE_Leave1_EUC = [CVE_Leave1_EUC,   err_per_case(:, 1)];
    CVE_Leave1_AFF = [CVE_Leave1_AFF,  err_per_case(:, 2)];
    CVE_Leave1_TPS3 = [CVE_Leave1_TPS3,  err_per_case(:, 3)];
    CVE_Leave1_TPS5 = [CVE_Leave1_TPS5,  err_per_case(:, 4)];
    CVE_Leave1_TPS7 = [CVE_Leave1_TPS7,  err_per_case(:, 5)];
    CVE_Leave1_Kernel = [CVE_Leave1_Kernel,  err_per_case(:, 6)];
end


CVE_Nfold_EUC = [];
CVE_Nfold_AFF = [];
CVE_Nfold_TPS3 = [];
CVE_Nfold_TPS5 = [];
CVE_Nfold_TPS7 = [];
CVE_Nfold_Kernel = [];

for ii = 1 : length(num_correspondences_to_remove)
    ijcv_case = IJCV2022_CVE_Nfold{ii};
    rss_case = RSS2022_CVE_Nfold{ii};  
    err_per_case = (ijcv_case - rss_case) ./ ijcv_case;
    CVE_Nfold_EUC = [CVE_Nfold_EUC,   err_per_case(:, 1)];
    CVE_Nfold_AFF = [CVE_Nfold_AFF,  err_per_case(:, 2)];
    CVE_Nfold_TPS3 = [CVE_Nfold_TPS3,  err_per_case(:, 3)];
    CVE_Nfold_TPS5 = [CVE_Nfold_TPS5,  err_per_case(:, 4)];
    CVE_Nfold_TPS7 = [CVE_Nfold_TPS7,  err_per_case(:, 5)];
    CVE_Nfold_Kernel = [CVE_Nfold_Kernel,  err_per_case(:, 6)];    
end




if use_percentage
    
    CVE_Leave1_EUC = CVE_Leave1_EUC * 100;
    CVE_Leave1_AFF = CVE_Leave1_AFF * 100;
    CVE_Leave1_TPS3 = CVE_Leave1_TPS3 * 100;
    CVE_Leave1_TPS5 = CVE_Leave1_TPS5 * 100;
    CVE_Leave1_TPS7 = CVE_Leave1_TPS7 * 100;
    CVE_Leave1_Kernel = CVE_Leave1_Kernel * 100;
    
    CVE_Nfold_EUC = CVE_Nfold_EUC * 100;
    CVE_Nfold_AFF = CVE_Nfold_AFF * 100;
    CVE_Nfold_TPS3 = CVE_Nfold_TPS3 * 100;
    CVE_Nfold_TPS5 = CVE_Nfold_TPS5 * 100;
    CVE_Nfold_TPS7 = CVE_Nfold_TPS7 * 100;
    CVE_Nfold_Kernel = CVE_Nfold_Kernel * 100;
    
end



colors = parula(6);


color1 = colors(1, :);
color2 = colors(2, :);
color3 = colors(3, :);
color4 = colors(4, :);
color5 = colors(5, :);


xvalues = num_correspondences_to_remove;
xxvalues = kron(num_correspondences_to_remove, ones(1, MC_Trials));

MethodNames ={"EUC", "AFF", "TPS(3)", "TPS(5)", "TPS(7)", "Kernel"};

label_font_size = 7;

close all;

hFig = figure("Name", datasetName, "Position", [200, 200, 500, 250]);
tfig = tiledlayout(1, 2, "TileSpacing", "tight");

colors = hsv(6);

ax1 = nexttile;    
ax1.FontSize = 5;
ax1.LineWidth = 0.5;
ax1.TickDir = 'in';

scatter_point_size = 8;

hold on;
h1 = errorbar(xvalues, mean(CVE_Leave1_AFF), std(CVE_Leave1_AFF), 'Color', color1);
c1 = scatter(xxvalues, CVE_Leave1_AFF(:), scatter_point_size, color1, '.');
h2 = errorbar(xvalues, mean(CVE_Leave1_TPS3), std(CVE_Leave1_TPS3), 'Color', color2);
c2 = scatter(xxvalues, CVE_Leave1_TPS3(:), scatter_point_size, color2, '.');
h3 = errorbar(xvalues, mean(CVE_Leave1_TPS5), std(CVE_Leave1_TPS5), 'Color', color3);
c3 = scatter(xxvalues, CVE_Leave1_TPS5(:), scatter_point_size, color3, '.');
h4 = errorbar(xvalues, mean(CVE_Leave1_TPS7), std(CVE_Leave1_TPS7), 'Color', color4);
c4 = scatter(xxvalues, CVE_Leave1_TPS7(:), scatter_point_size, color4, '.');
h5 = errorbar(xvalues, mean(CVE_Leave1_Kernel), std(CVE_Leave1_Kernel), 'Color', color5);
c5 = scatter(xxvalues, CVE_Leave1_Kernel(:), scatter_point_size, color5, '.');
yline(ax1, 0, 'r-', {'Above this line', 'the proposed is better'}, 'FontSize', 6, 'LabelHorizontalAlignment', 'right', 'interpreter', 'latex');
hold off;
title("CVE - leave 1", "FontSize", label_font_size, 'FontWeight', 'normal');
xticks(xvalues);
if use_percentage
    ytickformat('percentage')
end
axis tight; box on;


ax2 = nexttile;    
ax2.FontSize = 5;
ax2.LineWidth = 0.5;
ax2.TickDir = 'in';

hold on;
hh1 = errorbar(xvalues, mean(CVE_Nfold_AFF), std(CVE_Nfold_AFF), 'Color', color1);
cc1 = scatter(xxvalues, CVE_Nfold_AFF(:), scatter_point_size, color1, '.');
hh2 = errorbar(xvalues, mean(CVE_Nfold_TPS3), std(CVE_Nfold_TPS3), 'Color', color2);
cc2 = scatter(xxvalues, CVE_Nfold_TPS3(:), scatter_point_size, color2, '.');
hh3 = errorbar(xvalues, mean(CVE_Nfold_TPS5), std(CVE_Nfold_TPS5), 'Color', color3);
cc3 = scatter(xxvalues, CVE_Nfold_TPS5(:), scatter_point_size, color3, '.');
hh4 = errorbar(xvalues, mean(CVE_Nfold_TPS7), std(CVE_Nfold_TPS7), 'Color', color4);
cc4 = scatter(xxvalues, CVE_Nfold_TPS7(:), scatter_point_size, color4, '.');
hh5 = errorbar(xvalues, mean(CVE_Nfold_Kernel), std(CVE_Nfold_Kernel), 'Color', color5);
cc5 = scatter(xxvalues, CVE_Nfold_Kernel(:), scatter_point_size, color5, '.');
yline(ax2, 0, 'r-', {'Above this line', 'the proposed is better'}, 'FontSize', 6, 'LabelHorizontalAlignment', 'right', 'interpreter', 'latex');
hold off;
title(sprintf("CVE - %d fold", Nfold), "FontSize", label_font_size, 'FontWeight', 'normal');
xticks(xvalues);
if use_percentage
    ytickformat('percentage')
end
axis tight; box on;

lgd = legend([hh1, hh2, hh3, hh4, hh5], MethodNames(2:6), 'FontSize', label_font_size, 'Orientation', 'Horizontal', 'box', 'off');
lgd.Layout.Tile = 'north';

xlabel(tfig, "number of missing correspondences", "FontSize", label_font_size);
exportgraphics(hFig, [fdir, 'Lambda_Methods_by_Missing_Correspondences_', datasetName, '.pdf'], 'ContentType','vector');







function [cve_rigid, cve_affine, cve_tps3, cve_tps5, cve_tps7, cve_kernel] = runAllMethods (DataCell, VisVecCell, ...
    TPS_smooth_param, TPS_control_points_dim, ...
    Kernel_p_quantile, Kernel_smooth_param, ...
    flagJointlyOptimizeLamdaRigidTransformation, Nfold)

    
%     GPA_RIGID = RigidGPA;
%     GPA_RIGID.addDataByArray(DataCell, VisVecCell);
%     GPA_RIGID.verbose = false;
%     GPA_RIGID.run();
%     cve_rigid = GPA_RIGID.run_N_Fold_CrossValidation(Nfold);
    cve_rigid = 0;


    
    GPA_AFFINE = DefGPA('AFFINE');
    GPA_AFFINE.verbose = false;
    GPA_AFFINE.dimFeatureSpace = 0;
    GPA_AFFINE.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
    GPA_AFFINE.addDataByArray(DataCell, VisVecCell);
    GPA_AFFINE.run (TPS_smooth_param, flagJointlyOptimizeLamdaRigidTransformation);
    GPA_AFFINE.rsdError();
    cve_affine = GPA_AFFINE.run_N_Fold_CrossValidation(Nfold);
    
    
    GPA_TPS3 = DefGPA('TPS');
    GPA_TPS3.verbose = false;
    GPA_TPS3.dimFeatureSpace =3^TPS_control_points_dim;
    GPA_TPS3.smoothParam = TPS_smooth_param;
    GPA_TPS3.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
    GPA_TPS3.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
    GPA_TPS3.addDataByArray(DataCell, VisVecCell);
    GPA_TPS3.run (TPS_smooth_param, flagJointlyOptimizeLamdaRigidTransformation);
    GPA_TPS3.rsdError();
    cve_tps3 = GPA_TPS3.run_N_Fold_CrossValidation(Nfold);

    
    GPA_TPS5 = DefGPA('TPS');
    GPA_TPS5.verbose = false;
    GPA_TPS5.dimFeatureSpace =5^TPS_control_points_dim;
    GPA_TPS5.smoothParam = TPS_smooth_param;
    GPA_TPS5.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
    GPA_TPS5.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
    GPA_TPS5.addDataByArray(DataCell, VisVecCell);
    GPA_TPS5.run (TPS_smooth_param, flagJointlyOptimizeLamdaRigidTransformation);
    GPA_TPS5.rsdError();
    cve_tps5 = GPA_TPS5.run_N_Fold_CrossValidation(Nfold);

    
    
    GPA_TPS7 = DefGPA('TPS');
    GPA_TPS7.verbose = false;
    GPA_TPS7.dimFeatureSpace = 7^TPS_control_points_dim;
    GPA_TPS7.smoothParam = TPS_smooth_param;
    GPA_TPS7.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
    GPA_TPS7.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
    GPA_TPS7.addDataByArray(DataCell, VisVecCell);
    GPA_TPS7.run (TPS_smooth_param, flagJointlyOptimizeLamdaRigidTransformation);
    GPA_TPS7.rsdError();
    cve_tps7 = GPA_TPS7.run_N_Fold_CrossValidation(Nfold);
    
    
    GPA_KERNEL = KernelGPA;
    GPA_KERNEL.verbose = false;
    GPA_KERNEL.smoothParam = Kernel_smooth_param;
    GPA_KERNEL.quantilePercentParam= Kernel_p_quantile;
    GPA_KERNEL.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
    GPA_KERNEL.centerDataAhead = false; % 1 -- center data ahead.  0 -- don't center data ahead
    GPA_KERNEL.addDataByArray(DataCell, VisVecCell);
    GPA_KERNEL.run ();
    GPA_KERNEL.rsdError();
    cve_kernel = GPA_KERNEL.run_N_Fold_CrossValidation(Nfold);

end





function [DataCell, VisVecCell, all_vis] = RemoveCorrespondences (DataCell, VisVecCell, N)

numPointClouds = length(DataCell);
numPoints = size(DataCell{1}, 2);

all_vis = zeros(1, numPoints);
for ii = 1 : length(VisVecCell)
    all_vis = all_vis + VisVecCell{ii};
end

for ii = 1 : N
    reps = 1;
    while (reps < 10000)
        rn = max(1, round(rand * numPointClouds * numPoints));
        cloud_idx = ceil(rn/numPoints);
        point_idx = rn - (cloud_idx-1)*numPoints;
        if (VisVecCell{cloud_idx}(point_idx) && all_vis(point_idx) > 2)
            VisVecCell{cloud_idx}(point_idx) = false;
            all_vis(point_idx) = all_vis(point_idx) - 1;
            reps = inf;
        end
        reps = reps + 1;
    end
end

% check results
all_vis = zeros(1, size(DataCell{1}, 2));
for ii = 1 : length(VisVecCell)
    all_vis = all_vis + VisVecCell{ii};
end
all_vis;

end





