function [Data] = cleanInvisiblePoints(Data)
% cleanInvisiblePoints: cleans data of invisible points by assigning zeros
% to point correspondences and GT
    viz = Data.v;
    nFiles = size(viz,1);
    nFeats = size(viz,2);
    
    for iFile = 1:nFiles
        p = Data.p(iFile).p;
        P = Data.Pgth(iFile).P;
        for iFeat = 1:nFeats
            if(viz(iFile, iFeat) == 0)
                p(:,iFeat) = 0;
                P(:,iFeat) = 0;
            end
        end
        Data.p(iFile).p = p;
        Data.Pgth(iFile).P = P;
    end

end