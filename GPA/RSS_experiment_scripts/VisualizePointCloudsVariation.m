function error_statistics = VisualizePointCloudsVariation (PointCloudsCell, VisVecCell)

hold on;

num_pointClouds = length(PointCloudsCell);

array_num_points = zeros(1, num_pointClouds);
for ii = 1 : num_pointClouds
    array_num_points(ii) = length(PointCloudsCell{ii});
end
if ( min(array_num_points) ~= max(array_num_points) )
    fprintf(2, 'Different number of points in each point-clouds \n');
else
    num_points = array_num_points(1);
end

dim = size(PointCloudsCell{1}, 1);
ref_point_cloud = zeros(dim, num_points);
sum_vis = zeros(1, num_points);






if (true)
    C = jet(num_points+5);
%    C = winter(num_points);
%    C = hsv(num_points+30);
%    C = colorcube(num_points+50);
    for ii = 1 : num_pointClouds
        Points = PointCloudsCell{ii};
        vis = logical(VisVecCell{ii});
        if (size(Points, 1) == 2)
            scatter(Points(1, vis), Points(2, vis), 1,  C(vis, :), 'filled', "Marker", 'o');
        end
        if (size(Points, 1) == 3)
            scatter3(Points(1, vis), Points(2, vis), Points(3, vis), 1, C(vis, :), 'filled', "Marker", 'o');
        end
        ref_point_cloud = ref_point_cloud + Points .* vis;
        sum_vis = sum_vis + vis;
    end
else
    C = cool(num_pointClouds+2);
    for ii = 1 : num_pointClouds
        color_point_cloud = C(ii, :);
        Points = PointCloudsCell{ii};
        vis = logical(VisVecCell{ii});
        if (size(Points, 1) == 2)
            scatter(Points(1, vis), Points(2, vis), 1,  color_point_cloud, "Marker", 'x');
        end
        if (size(Points, 1) == 3)
            scatter3(Points(1, vis), Points(2, vis), Points(3, vis), 1, color_point_cloud, "Marker", 'x');
        end
        ref_point_cloud = ref_point_cloud + Points .* vis;
        sum_vis = sum_vis + vis;
    end
end


for jj = 1 : length(sum_vis)
    if (sum_vis(jj) > 0)
        ref_point_cloud(:, jj) = ref_point_cloud(:, jj) ./ sum_vis(jj);
    end
end

squared_accum_matr = 0 * ref_point_cloud;
for ii = 1 : num_pointClouds
    Points = PointCloudsCell{ii};
    vis = logical(VisVecCell{ii});
    tmp = (Points - ref_point_cloud) .* vis;
    squared_accum_matr = squared_accum_matr + tmp .* tmp;    
end
points_std_dev = sqrt( sum(squared_accum_matr, 1) ./ sum_vis );


color_ref_edge = '#18C1FF';
color_ref_face  = 'none';

if (size(ref_point_cloud, 1) == 2)
    scatter(ref_point_cloud(1, :), ref_point_cloud(2, :), points_std_dev,  'MarkerEdgeColor', color_ref_edge, 'MarkerFaceColor', color_ref_face);
end
if (size(ref_point_cloud, 1) == 3)
    scatter3(ref_point_cloud(1, :), ref_point_cloud(2, :), ref_point_cloud(3, :), points_std_dev, 'MarkerEdgeColor', color_ref_edge, 'MarkerFaceColor', color_ref_face);
end


error_statistics = sqrt( sum(sum(squared_accum_matr, 1)) / sum(sum_vis) );


if (dim == 3)
    view(3);
end

hold off;

end

