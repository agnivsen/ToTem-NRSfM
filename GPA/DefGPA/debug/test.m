% clear all
% close all
% clc
% 
% 
% genNewData = true;
% 
% inputFile = 'dataset/DeformedMesh/dmesh0.txt';
% dirPath = 'dataset/DeformedMeshIdentical';
% 
% if (genNewData)
%     delete ([dirPath,'/*']);
%     mkdir (dirPath);
%     for cnt = 1 : 6
%         axis = randn(3,1);
%         axis = axis / norm(axis);
%         Rot = SO3.Exp(0.5 * cnt * axis);
%         Rot = 10* rand(3,3);
%         genData (inputFile,  [dirPath, '/rot', num2str(cnt), '.txt'], Rot);
%     end
% end
% 
% 
% WP = WarpProcrustesAnalysis;
% 
% WP.dimFeatureSpace = 20;
% 
% WP.readDir (dirPath);
% 
% WP.smoothParam = 1   % * WP.numPoints
% 
% WP.run ();
% 
% WP.writeShape ([dirPath, '/optReference']);





clear all
close all
clc


genNewData = true;

dirPath = 'dataset/Karim_v2/DeformedMesh';

%dirPath = 'dataset/DeformedMeshIdentical';
% if (genNewData)
%     delete ([dirPath,'/*']);
%     mkdir (dirPath);
%     for cnt = 1 : 6
%         axis = randn(3,1);
%         axis = axis / norm(axis);
%         Rot = SO3.Exp(0.5 * cnt * axis);
%         % Rot = 10* rand(3,3);
%         genData (inputFile,  [dirPath, '/rot', num2str(cnt), '.txt'], Rot);
%     end
% end


WP = DGPA_Warp;

WP.dimFeatureSpace = 20;

WP.readDir (dirPath);

WP.smoothParam = 1   % * WP.numPoints

WP.run ();

WP.writeShape ([dirPath, '/optReference']);





function genData (InputFile, OutputFile, Rotation)


dataMatr = [];
fid = fopen (InputFile, 'r');
lstr = fgetl (fid);
pcnt = 0;

while ischar (lstr)
    
    if (size (lstr, 2) > 0)
        
        pcnt = pcnt + 1;
        
        dataMatr(:, pcnt) = str2double(split(lstr));
        
    end
    
    lstr = fgetl (fid);
    
end

fclose (fid);

% write data

dataMatr = Rotation * dataMatr;

dataMatr = dataMatr + 5 * rand (size(dataMatr, 1), size(dataMatr,2));


fwt = fopen (OutputFile, 'w');
formatSpec = '%f %f %f\n';

for ii = 1 : size(dataMatr, 2)
    fprintf(fwt, formatSpec, dataMatr(:, ii));
end

fclose(fwt);

end


