function [eSquaredSummed, eListSquared, eListSigned] = evaluateTotalCostCombined(directionVectors, nng, depth, geodesic, W, visibilityMatrix)

    nFiles =  size(directionVectors,2);
    
    nPts = size(directionVectors{1},2);
    
    nNeighbors = size(nng, 2);
    
    cost = CostNRSfM(W);
    
    listE = []; 
    
    for iFile = 1:nFiles
        for iAlpha = 1:nPts
            for iBeta = 1:nNeighbors
                if(nng(iAlpha, iBeta) > 0)
                    
                    iNeigh = nng(iAlpha, iBeta);
                    
                    if((visibilityMatrix(iFile, iAlpha) == 1) && (visibilityMatrix(iFile, iNeigh) == 1))
                        alpha = depth{iFile}(1, iAlpha);
                        beta = depth{iFile}(1, iNeigh);

                        qJ = directionVectors{iFile}(:,iAlpha);
                        qQ = directionVectors{iFile}(:,iNeigh);
                        [theta] = angleBetweenVectors(qJ,qQ);             

                        d = geodesic{iAlpha, iNeigh};

                        if( isempty(d) || d <= 0 )
                            d = geodesic{iNeigh, iAlpha}; % since geodesics  A -> B  = B -> A
                            if(isempty(d) || d <= 0)
                                fprintf('evaluateTotalCostCombined: geodesic distance between %d and %d-th keypoints in %d-th file is %d\n', iAlpha, iNeigh, iFile, d);
                                error('evaluateTotalCostCombined: invalid geodesic distance received (either negative distance or empty cell)!');
                            end
                        end


                         [C] =  cost.costSigned(alpha, beta, theta, d);   

                         listE = [listE C]; 
                    end                    
                end
            end
        end
    end
    
    eSquaredSummed = sum(listE.^2)/numel(listE);
    listE = listE./sqrt(numel(listE));
    eListSigned = sum(listE);
    eListSquared = eListSigned.^2;

end