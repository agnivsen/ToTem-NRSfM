function hFig = VisualizePointClouds(GPA_Handle)

nsize = GPA_Handle.numPointClouds;

hFig = figure;

tfig = tiledlayout(2, nsize+1);
tfig.TileSpacing = 'compact';
% tfig.TileSpacing = 'tight';
% tfig.TileSpacing = 'none';


% ----- the original point-clouds ----

for ii = 1 : nsize
    
    data = GPA_Handle.PointClouds(ii).Data;
    vis = GPA_Handle.PointClouds(ii).Vis;
    trans = GPA_Handle.PointClouds(ii).TransPrior;
    data_vis = data(:, logical(vis)) - trans;
    
    ax = nexttile;
    if (size(data, 1) == 2)
        scatter(data_vis(1,:), data_vis(2,:));
    end
    if (size(data, 1) == 3)
        scatter3(data_vis(1,:), data_vis(2,:), data_vis(3,:));
    end
    title(ax, ['data ', num2str(ii)]);
    
end

ax = nexttile;
title(ax, 'reference');


% ----- the transformed point-clouds ----

for ii = 1 : nsize
    
    data = GPA_Handle.PointClouds(ii).Data;
    vis = GPA_Handle.PointClouds(ii).Vis;
    trans = GPA_Handle.PointClouds(ii).TransPrior;
    data_vis = data(:, logical(vis)) - trans;
    
    data_vis = GPA_Handle.transformPoints(data_vis, ii);
    
    ax = nexttile;
    if (size(data, 1) == 2)
        scatter(data_vis(1,:), data_vis(2,:));
    end
    if (size(data, 1) == 3)
        scatter3(data_vis(1,:), data_vis(2,:), data_vis(3,:));
    end
    
end

ax = nexttile;
ref_shape = GPA_Handle.mShape;
if (size(data, 1) == 2)
    scatter(ref_shape(1,:), ref_shape(2,:));
end
if (size(data, 1) == 3)
    scatter3(ref_shape(1,:), ref_shape(2,:), ref_shape(3,:));
end



set(0,'units','pixels');
Pix_SS = get(0,'screensize');
hFig.Position = [1, Pix_SS(4), Pix_SS(3), Pix_SS(4)/3];

end

