clear all
close all
clc

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




fdir = '../../dataset/TOPACS/landmarks/';

num_frames = 6;
for i = 1 : num_frames
    PointsCell{i} = readmatrix([fdir,'points', num2str(i-1), '.csv']);
end
pairs = readmatrix([fdir,'pairs.csv']);


if(1)
    pos_arr = { [0, 800, 400, 500], ...
        [600, 800, 400, 500], ...
        [1200, 800, 400, 500], ...
        [0, 0, 400, 500], ...
        [600, 0, 400, 500], ...
        [1200, 0, 400, 500] };
    
    for ii = 1 : num_frames
        hFig = figure(ii);
        pts = PointsCell{ii}(:, 1:3);
        scatter3(pts(:,1), pts(:,2), pts(:,3), '.');
        hFig.Position = pos_arr{ii};
        axis equal;
    end
end




return;

[DataCell, VisVecCell] = ExtractTOPACS(fdir);


% find points visible in all shapes
covis = ones(1, 20000);
for kk = 1 : num_frames
    covis = covis .* VisVecCell{kk};
end

covis = logical(covis);
    
for kk = 1 : num_frames
    DataCell{kk} = DataCell{kk}(:, covis);
    VisVecCell{kk} = VisVecCell{kk}(covis);
end


nsize = size(DataCell{1},2)




if (1)

EUC_ALL = GPA_Functor('EUC-ALL');
EUC_ALL.verbose = 'on';
EUC_ALL.addDataByArray(DataCell, VisVecCell);
EUC_ALL.run();
EUC_ALL.rsdError();
EUC_ALL



AFF_ALL = GPA_Functor('AFF-ALL');
AFF_ALL.verbose = 'on';
AFF_ALL.addDataByArray(DataCell, VisVecCell);
AFF_ALL.run();
AFF_ALL.rsdError();
AFF_ALL




GPA_AFFINE = DefGPA('AFFINE');
GPA_AFFINE.addDataByArray(DataCell, VisVecCell);
GPA_AFFINE.run ()
GPA_AFFINE.rsdError();
GPA_AFFINE




GPA_TPS_case3 = DefGPA('TPS');
GPA_TPS_case3.addDataByArray(DataCell, VisVecCell);
GPA_TPS_case3.dimFeatureSpace = 5^3;
GPA_TPS_case3.run (0.1);
GPA_TPS_case3.rsdError();
GPA_TPS_case3



end















function [points, identifiers] = ReadLiTSFCSV (fileName)

points = [];

identifiers = '';

fid = fopen (fileName);

lstr = fgetl (fid);

while ischar (lstr)
    
    if (size (lstr, 2) > 0) % skip empty lines
        
        strArr = split (lstr, ',');
        leadingChar = strArr{1};
        
        if ~isempty(str2num(leadingChar))
            
            x = str2double(strArr{2});
            y = str2double(strArr{3});
            z = str2double(strArr{4});
            
            id_str = strArr{12};
            
            points = [points,  [x; y; z;] ];
            
            identifiers = [identifiers, ' ', id_str];
            
        end

    end
    
    lstr = fgetl (fid);
    
end

fclose (fid);

end









function [DataCell, VisVecCell] = ExtractTOPACS (fdir)

num_frames = 6;

for i = 1 : num_frames
    PointsCell{i} = readmatrix([fdir,'points', num2str(i-1), '.csv']);
end
pairs = readmatrix([fdir,'pairs.csv']);


% construct adjacent arrays

for ii = 1 : num_frames
    AdjacentArray{ii} = zeros(20000, num_frames);
end

for kk = 1 : size(pairs, 1)
    
    val = pairs(kk, :);
    
    imgID1 = val(1)+1;
    ptsID1 = val(2)+1;
    imgID2 = val(3)+1;
    ptsID2 = val(4)+1;
    distance = val(5);
    probability = val(6);
    
    
    if (probability < 0.94 || distance > 10)
        continue;
    end
    
    
    fprintf('val = [%d, %d, %d, %d, %f, %f]\n', val(1), val(2), val(3), val(4), val(5), val(6));
    % construact adjacent list
    
    v1 = AdjacentArray{imgID1}(ptsID1, :);
    v2 = AdjacentArray{imgID2}(ptsID2, :);

    
%     v1(imgID1) = ptsID1;
%     v2(imgID2) = ptsID2;

    % deal with two hypothesis here?
    if (v1(imgID2) ==0)
        v1(imgID2) = ptsID2;
    elseif (v1(imgID2) ~= ptsID2)
        fprintf(2, 'ambiguous pairs\n');
        v1(imgID2) = 0;
    end
    if (v2(imgID1) == 0)
        v2(imgID1) = ptsID1;
    elseif (v2(imgID1) ~= ptsID1)
        fprintf(2, 'ambiguous pairs\n');
        v2(imgID1)  = 0;
    end
    
    
    
    v1(imgID2) = ptsID2;
    v2(imgID1) = ptsID1;

    fprintf('[%d] v1 = [%d %d %d %d %d %d],  [%d] \n', imgID1, v1, ptsID1);
    fprintf('[%d] v2 = [%d %d %d %d %d %d],  [%d] \n\n', imgID2, v2, ptsID2);
    
    
    AdjacentArray{imgID1}(ptsID1, :) = v1;
    AdjacentArray{imgID2}(ptsID2, :) = v2;
    
end



% initialize storage

DataCell = cell(1, num_frames);
VisVecCell = cell(1, num_frames);
for kk = 1 : num_frames
    DataCell{kk} = zeros(3, 20000);
    VisVecCell{kk} = zeros(1, 20000);
end



for ii = 1 : 20000
    
    matx = zeros(num_frames);
    
    vtmp = AdjacentArray{1}(ii, :);
    
    for loop = 1 : 6
        for kk = 1 : num_frames
            if (vtmp(kk) >0)
                matx(kk, :) = AdjacentArray{kk}(vtmp(kk), :);
            end
        end
        % sychronize adjacent vectors
        vnew = max(matx);
        if (vnew == vtmp)
            break;
        else
            vtmp = vnew;
        end
    end
    

    vtmp = max (matx);
    
    if (nnz(vtmp) == 6)
        matx
    end
    
    
    for kk = 1 : num_frames
        if (vtmp(kk) > 0)
            DataCell{kk}(:, ii) = PointsCell{kk}(vtmp(kk), [1,2,3])';
            VisVecCell{kk}(ii) = 1;
        end
    end
    
end

end






    