function error_statistics = VisualizePointCloudsVariation (TestPointsDataCell, TestPointsVisVecCell)


[points_std_dev, mean_template, RMSE_VE] = ComputePointwiseConsistency (TestPointsDataCell, TestPointsVisVecCell);

dim = size(mean_template, 1);

num_points = length(points_std_dev);

C = colorcube(num_points+5);

hold on;

for ii = 1 : length(TestPointsDataCell)
    vis = logical(TestPointsVisVecCell{ii});
    Points = TestPointsDataCell{ii} .* vis;
    if (dim == 2)
        scatter(Points(1, vis), Points(2, vis), 1,  C(vis, :), 'filled', "Marker", 'o');
    end
    if (dim == 3)
        scatter3(Points(1, vis), Points(2, vis), Points(3, vis), 1, C(vis, :), 'filled', "Marker", 'o');
    end
end



%color_ref_edge = '#18C1FF';
color_ref_edge = '#111111';
color_ref_face  = [0.7, 0.7, 0.7];

if (dim == 2)
    scatter(mean_template(1, :), mean_template(2, :), points_std_dev, points_std_dev, 'filled'); %, 'o', "MarkerEdgeColor", C);
end
if (dim == 3)
    scatter3(mean_template(1, :), mean_template(2, :), mean_template(3, :), points_std_dev, points_std_dev, 'filled');% , 'o', "MarkerEdgeColor", C);
end
colormap(flipud(gray));
%colormap (0.95*flipud(pink));
caxis([0, 60]);
%cb = colorbar('FontSize',10, 'Location', 'northoutside');


error_statistics.min = min(points_std_dev);
error_statistics.max = max(points_std_dev);
error_statistics.mean = mean(points_std_dev);
error_statistics.RMSE_ve = RMSE_VE;

if (dim == 3)
    view(3);
end

hold off;

end

