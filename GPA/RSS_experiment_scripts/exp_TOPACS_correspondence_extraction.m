clear all
close all
clc


fdir = '../../dataset/TOPACS/landmarks/';


distance_tol = 8;
probability_tol = 0.94;
corresp_minimum_shapes = 2;

[DataCell, VisVecCell] = ProcessTOPACS (fdir, distance_tol, probability_tol, corresp_minimum_shapes)

 

 % count the number of visibilities for each correspondece
 num_point_clouds = length(VisVecCell);
 num_points = length(VisVecCell{1});
 VisiblityMatrix = zeros(num_point_clouds,  num_points);
for ii = 1 : num_point_clouds
    vis = VisVecCell{ii};
    VisiblityMatrix(ii, :) = reshape(vis, 1, num_points);
end
numVis = sum(VisiblityMatrix, 1);
fprintf(1, '\n');
for cnt_vis = 2 : 6
    fprintf(1, '%d points occur in exactly %d point-clouds\n', sum(numVis == cnt_vis), cnt_vis);
end


% There are 6 volumes containing 20,000 key points, and approximately 300,000 pairs
% 
% * pointsXX.csv.gz: list of key points for each image. Each line of the file corresponds to a key point, and contains 54 values:
% 
% - The 3D coordinates of the point
% - the point scale
% - the sign of the Laplacian
% - the value of the detector
% - the SURF3D descriptor vector (48 values)
% 
% * pairs.csv.gz: list of pairs. Each row contains 6 values ​​for each pair:
% 
% - the id of image 1
% - the id of point 1 in image 1
% - the Id of image 2
% - the Id of point 2 in image 2
% - the distance between the points after registration
% - the probability that the pair is a true pair according to FROG. This probability is calculated from the statistics in Figure 1.
% 
% Pairs are sorted by distance (5th value in the row)
% 
% BE CAREFUL, these are in fact half-pairs, so each pair appears twice in the list, because the probability has two possible values, one for each image. If you also need the coordinates of the readjusted points, I can add them, up to you.
% 
% I can send you the original images too if needed. I have attached to this message a capture of the average image after registration by FROG. We can see that one of the images is smaller than the others (upper body only), and that some patients have their arms raised, others not. In short, very real data!





function [DataCell, VisVecCell] = ProcessTOPACS (fdir, distance_tol, probability_tol, corresp_minimum_shapes)

pairs = SelectParisFromTOPACS (fdir, distance_tol, probability_tol);

NodeOffsetEncoding = 20000;

% a cell, with each elements a 2 rows matrix:   first-row for ImgeIDs;  second-row for PointIDs
bins = AssociatePairsGlobally (pairs, NodeOffsetEncoding);

clear pairs;


% read all point-clouds
numPointClouds = 6;
for ii = 1 : numPointClouds
        tmp = readmatrix([fdir,'points', num2str(ii-1), '.csv']);
        PointCloudsCell{ii} = tmp(:, 1:3)';
end
clear tmp;


% remove ambigous associations
assMatrix = [];
for ii = 1 : length(bins)
    
    ass_ambigous = bins{ii};
    ImgID = ass_ambigous(1, :);
    PtsID = ass_ambigous(2, :);

    ass = - ones(1, numPointClouds);
    
    for jj = 1 : numPointClouds
        pointIDs_in_image_jj = PtsID(ImgID == jj);
        if (size(pointIDs_in_image_jj, 2) > 0)
            ass(jj) = pointIDs_in_image_jj(1);
            % in case there are multiple ambigous associations
            if (size(pointIDs_in_image_jj, 2) > 1)
                % use the central point to which the overall Euclidean distance is smallest
                points = PointCloudsCell{jj}(:, pointIDs_in_image_jj);   % ambigous points
                idx = GetCentralPoint(points);
                ass(jj) = pointIDs_in_image_jj(idx);
            end
        end
    end
    
    % make sure the correspondence occurs at least in two point-clouds
    if (nnz(ass > 0) > corresp_minimum_shapes)
        assMatrix = [assMatrix;  ass];
    end
    
end


% save processed association marix
writematrix(assMatrix, [fdir, 'CorrespondenceMatrix', '.txt'], "Delimiter", " ");
fprintf(2, ['save correspondence association matrix to path:\n',  fdir, 'CorrespondenceMatrix', '.txt',  '\n']);

% construct final ouput, from pointclouds and association-matrix
num_corrsp = size(assMatrix, 1);

for ii = 1 : numPointClouds
    v = assMatrix(:, ii);
    vis = logical(v > 0);    
    data = zeros(3, num_corrsp);   
    data(:, vis) = PointCloudsCell{ii}(:, v(vis));
    DataCell{ii} = data;
    VisVecCell{ii} = vis';
end

end




function idx = GetCentralPoint(points)
    nsize = size(points, 2);
    M = zeros(nsize);
    for ii = 1 : nsize
        for jj = (ii+1) : nsize
            dp = points(:, ii) - points(:, jj);
            M(ii, jj) = dp' * dp;
        end
    end
    M = M + M';
    squared_distance = sum(M);
    [~, idx] = min(squared_distance);
end




function pairs = SelectParisFromTOPACS (fdir, distance_tol, probability_tol)

pairs = readmatrix([fdir,'pairs.csv']);

idx = [];

for kk = 1 : 2 : size(pairs, 1)
    
    val_forward = pairs(kk, :);
    distance_forward = val_forward(5);
    probability_forward = val_forward(6);
    
    val_backward = pairs(kk+1, :);
    distance_backward = val_backward(5);
    probability_backward = val_backward(6);
    
    distance = max(distance_forward, distance_backward);
    probability = min(probability_forward, probability_backward);
    
    if (distance < distance_tol && probability > probability_tol)
        idx = [idx; kk];
    end
    
end

pairs = pairs(idx, :);

end





function bins = AssociatePairsGlobally (pairs, NodeOffsetEncoding)

    nsize = length(pairs);
    
    % Node ID: encode point id and image id together
    edges = zeros(nsize, 2);
    for ii = 1 : nsize
        val = pairs(ii, :);
        imgID1 = val(1)+1; ptsID1 = val(2)+1;
        imgID2 = val(3)+1; ptsID2 = val(4)+1;
        node_id1 = NodeOffsetEncoding * (imgID1-1) + ptsID1;
        node_id2 = NodeOffsetEncoding * (imgID2-1) + ptsID2;
        edges(ii, :) = [node_id1, node_id2];
    end
    
    % construct graph
    G = graph(edges(:, 1), edges(:, 2));
    
    % compute connected components
    bins = G.conncomp('OutputForm', 'cell');
    
    % select connect component with size greater than 3
    comp_idx_used = [];
    for ii = 1 : length(bins)
        vec = bins{ii};
        if (length(vec) > 3)
            comp_idx_used = [comp_idx_used, ii];
        end
    end
    bins = bins(comp_idx_used);

    % map node index back to image-id and point-id
    for ii = 1 : length(bins)
        vec = bins{ii};
        PtIDs = mod(vec, NodeOffsetEncoding);
        ImgIDs = (vec - PtIDs) ./ NodeOffsetEncoding;
        ImgIDs = ImgIDs + 1;
        bins{ii} = [ImgIDs;  PtIDs];
    end
    
end




    