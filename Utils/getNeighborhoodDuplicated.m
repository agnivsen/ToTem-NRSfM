function [neighborhood] = getNeighborhoodDuplicated(Data, nK)
    %% getNeighborhood compute neighborhood from Data, doesn't remove duplicated edges
    % PARAMETERS:
    %    Data: 'struct' type variable of feature correspondences in NRSfM
    %                   format
    %    nK: no. of neighbors required
    %
    nFeats = size(Data.p(1).p,2);

    assert(nK > 0,'ERROR: (getNeighborhood): nK should be greater than ZERO');
    assert(nK < nFeats,'ERROR: (getNeighborhood): nK should be less than total no. of feature correspondences');

    REFERENCE_FRAME_NNG = 1; % NNG is always computed from first frame

    neighborhood =  kNearestNeighbors(Data.p(REFERENCE_FRAME_NNG).p.', Data.p(REFERENCE_FRAME_NNG).p.', (nK+1));
    neighborhood(:,1) = [];

    return;

end