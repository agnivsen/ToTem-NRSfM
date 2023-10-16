function [directionVectors, depth] = getDirectionVectors(Data, K, depthInit, depthFromGT)
%% getDirectionVectors: obtain direction vectors for feature correspondences.
% PARAMETERS:
%    Data: 'struct' type variable of feature correspondences in NRSfM
%                   format
%    depthFromGT: if set to 1, depths are initialised from GT
%    K: [3x3] intrinsics
%
    nFiles = size(Data.p,2);
    nFeats = size(Data.p(1).p,2);
    
    Kinv = K\eye(3); % inverting intrinsics
    
    for iFiles = 1:nFiles
        if(size(Data.p(iFiles).p, 1) == 2)
            directionVectorsUnnormalised = Kinv*[Data.p(iFiles).p; ones(1, nFeats)]; % normalised img. pt.s
        else
            directionVectorsUnnormalised = Kinv*[Data.p(iFiles).p]; % normalised img. pt.s
        end
        directionVectors{iFiles} = normc(directionVectorsUnnormalised) .* Data.v(iFiles,:); % unit vector in R^3
        if(depthFromGT == 1)
            depth{iFiles} = VecNorm(Data.Pgth(iFiles).P);
        else
            depth{iFiles} = depthInit.*ones(1,size(directionVectorsUnnormalised,2)); % initialised to 'depth'
        end
    end
    
    return;
end