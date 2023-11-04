function [DataCell, VisVecCell, SampleCell] = readDataset (datasetName, ffpath)

if ~exist('ffpath', 'var')
    ffpath = '../dataset';
end

if strcmp(datasetName, 'Face')
    fdir =  [ffpath, '/face2D/'];
    load ([ffpath, '/face2D/ShapeDataFace.mat']);
    DataCell = D;
    n=size(DataCell,2);
    for i = 1 : n
        VisVecCell{i} = diag(V{i})';
        SampleCell{i} = [fdir, framelist{i}];
    end
end

if strcmp(datasetName, 'HandBag')
    fdir = [ffpath, '/handBag/imageData/'];
    load([ffpath, '/handBag/correspData2D/GT/Corresp.mat']);
    DataCell = Corresp.ptsImg;
    n=size(DataCell,2);
    m=size(DataCell{1},2);
    for i = 1 : n
        VisVecCell{i} = ones(1, m);
        DataCell{i} = DataCell{i}([1,2],:);
        SampleCell{i} = [fdir, Corresp.imNames{i}, '.bmp'];
    end
end

if strcmp(datasetName, 'PillowCover')
    fdir = [ffpath, '/pillowCover/imageData/'];
    load([ffpath, '/pillowCover/correspData2D/GT/Corresp.mat']);
    DataCell = Corresp.ptsImg;
    n=size(DataCell,2);
    m=size(DataCell{1},2);
    for i = 1 : n
        VisVecCell{i} = ones(1, m);
        DataCell{i} = DataCell{i}([1,2],:);
        SampleCell{i} = [fdir, Corresp.imNames{i}, '.bmp'];
    end
end

if strcmp(datasetName, 'DePoLL')
load([ffpath, '/DePoLL/DePoLL_markers_final.mat'], 'DataCell');
for i = 1 : length(DataCell)
    DataCell{i} = DataCell{i}(:, 1:60);
    VisVecCell{i} = ones(1, size(DataCell{i}, 2));
    SampleCell{i} = [ffpath, '/DePoLL/DePoLL_reference.png'];
end
DataCell = DataCell(7:13);
VisVecCell = VisVecCell(7:13);
end

if strcmp(datasetName, 'Liver')
    fdir =  [ffpath, '/SemisyntheticLiver/PointClouds/'];
    fhandles = dir([fdir, '/*.txt']);
    for i = 1 : size(fhandles, 1)
        fname = fhandles(i).name;
        D = readmatrix([fdir, fname]);
        DataCell{i} = D';
        VisVecCell{i} = ones(1, size(DataCell{i},2));
        SampleCell{i} = [ffpath, '/SemisyntheticLiver/snapshot0.png'];
    end
end


if strcmp(datasetName, 'StereoToyRug')
    fdir = [ffpath, '/stereoToyRug/'];
    X = load([fdir, 'X.txt']);
    for i = 1:651
        Q = X(3*i-2:3*i,:);
        DataCell{i} = 100 * Q;
        VisVecCell{i} = ones(1, size(DataCell{i}, 2));
        SampleCell{i}{1} = sprintf([fdir, 'I_left_%03d.jpg'],i);
        SampleCell{i}{2} = sprintf([fdir, 'I_right_%03d.jpg'],i);
    end
    DataCell = DataCell(1: 3:600);
    VisVecCell = VisVecCell (1:3:600);
    SampleCell = SampleCell(1:3:600);
end


if strcmp(datasetName, 'TOPACS')
    fdir = ['/home/agniva/Documents/papers/topacs6_v1p0/'];
    numPointClouds = 6;
    for ii = 1 : numPointClouds
        tmp = dlmread([fdir,'points', num2str(ii-1), '.csv']);
        PointCloudsCell{ii} = tmp(:, 1:3)';
    end
    assMatrix = dlmread([fdir, 'CorrespondenceMatrix', '.txt']);
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


if strcmp(datasetName, 'LiTS')
    fdir = [ffpath, '/LiTS/landmarks/'];
    fileList = dir([fdir, '*.fcsv']);
    num_files = length(fileList);
    DataCell = cell(0);
    IdentfierCell = cell(0);
    for k = 1 : num_files
        fileName = [fdir, fileList(k).name];
        [points, identifiers] = ReadLiTSFCSV (fileName);
        DataCell = [DataCell, points];
        IdentfierCell = [IdentfierCell, identifiers];
    end
    % check identifiers
    for k = 2 : num_files
        if (~strcmp (IdentfierCell{k-1}, IdentfierCell{k-1}))
            fprintf(2, 'correspondence error \n');
            fprintf('%s \n', IdentfierCell{k-1});
            fprintf('%s \n', IdentfierCell{k});
        end
    end
    for k = 1 : num_files
        VisVecCell{k} = ones(1, length(DataCell{k}));
    end
    SampleCell{1} = [fdir, 'full_body.png'];
    SampleCell{4} = [fdir, 'python_script.png']; 
    SampleCell{2} = [fdir, 'full_body_0.png'];
    SampleCell{3} = [fdir, 'full_body_2.png'];
end

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


