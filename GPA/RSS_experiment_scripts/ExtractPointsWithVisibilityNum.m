function [DataCell, VisVecCell] = ExtractPointsWithVisibilityNum (DataCell, VisVecCell, Narray)

num_point_clouds = length(VisVecCell);

array_num_points = zeros(1, num_point_clouds);
for ii = 1 : num_point_clouds
    array_num_points(ii) = length(VisVecCell{ii});
end
if ( min(array_num_points) ~= max(array_num_points) )
    fprintf(2, 'Different number of points in each point-clouds \n');
end
for ii = 1 : num_point_clouds
    array_num_points(ii) = length(DataCell{ii});
end
if ( min(array_num_points) ~= max(array_num_points) )
    fprintf(2, 'Different number of points in each point-clouds \n');
end

num_points = array_num_points(1);


VisiblityMatrix = zeros(num_point_clouds,  num_points);

for ii = 1 : num_point_clouds
    vis = VisVecCell{ii};
    VisiblityMatrix(ii, :) = reshape(vis, 1, num_points);
end


numVis = sum(VisiblityMatrix, 1);

selction_by_vis = zeros(1, num_points);

for ii = 1 : length(Narray)
    vis = (numVis == Narray(ii));
    selction_by_vis = selction_by_vis + vis;
end

selction_by_vis = logical(selction_by_vis);

for ii = 1 : num_point_clouds
    data = DataCell{ii}(:, selction_by_vis);
    DataCell{ii} = data;
    vis = VisVecCell{ii}(selction_by_vis);
    VisVecCell{ii} = vis;
end


end