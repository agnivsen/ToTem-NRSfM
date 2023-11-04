% this function computes the mean template and point-wise consistencies of transformed test points
function [points_std_dev, mean_template, RMSE_VE] = ComputePointwiseConsistency (TestPointsDataCell, TestPointsVisVecCell)

num_pointClouds = length(TestPointsDataCell);

array_num_points = zeros(1, num_pointClouds);
for ii = 1 : num_pointClouds
    array_num_points(ii) = length(TestPointsDataCell{ii});
end
if ( min(array_num_points) ~= max(array_num_points) )
    fprintf(2, 'Different number of points in each point-clouds \n');
else
    num_points = array_num_points(1);
end

dim = size(TestPointsDataCell{1}, 1);
mean_template = zeros(dim, num_points);
sum_vis = zeros(1, num_points);


for ii = 1 : num_pointClouds
    Points = TestPointsDataCell{ii};
    vis = logical(TestPointsVisVecCell{ii});
    mean_template = mean_template + Points .* vis;
    sum_vis = sum_vis + vis;
end

% construct the mean template
for jj = 1 : length(sum_vis)
    if (sum_vis(jj) > 0)
        mean_template(:, jj) = mean_template(:, jj) ./ sum_vis(jj);
    else
        mean_template(:, jj) = 0;
    end
end

% compute point-wise standard deviation
squared_accum_matr = 0 * mean_template;
for ii = 1 : num_pointClouds
    Points = TestPointsDataCell{ii};
    vis = logical(TestPointsVisVecCell{ii});
    tmp = (Points - mean_template) .* vis;
    squared_accum_matr = squared_accum_matr + tmp .* tmp;    
end
points_std_dev = sqrt( sum(squared_accum_matr, 1) ./ sum_vis );

RMSE_VE = sqrt ( sum(sum(squared_accum_matr, 1)) / sum(sum_vis) );

end

