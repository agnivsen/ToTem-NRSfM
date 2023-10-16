function [angleList] = precomputeAngles(directionVectors, visibilityMatrix, nng)
% precomputeAngles: computes angles between direction vectors, since these
%                                   angles are constant for a given data
%
% *RETURNS*:
% * _angleList_: a set of 'i' matrices, i being no. of
%                                           images/files. These matrices have the same structure as NNG
%                                           Contains the angles between neighbors, in RADIANS

    nPts = size(directionVectors{1},2);
    nFiles =  size(directionVectors,2);
    nNieghbors = size(nng, 2);

    for iFile = 1:nFiles
        angles = zeros(size(nng));
        for iPt = 1:nPts
            if(visibilityMatrix(iFile, iPt) == 1)
                vecJ = directionVectors{iFile}(:,iPt);
                for iNeighbor = 1:nNieghbors
                    indexQ = nng(iPt, iNeighbor);
                    if(visibilityMatrix(iFile, indexQ) == 1)
                        vecQ = directionVectors{iFile}(:,indexQ);
                        theta = angleBetweenVectors(vecJ, vecQ);
                        angles(iPt, iNeighbor) = theta;
                    end
                end
            end
        end
        angleList{iFile} = angles;
    end
end