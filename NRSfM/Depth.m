classdef Depth <  handle

    properties(Access = private)
            directionVectors = [];
            nng = [];
            depth = [];
            geodesic = [] ;
            nNeighbors = 0;
            isDataSet = false;
            currAlpha = 0;
            currFile = 0;
            MAX_DEPTH = 2.0;
            MIN_DEPTH = 0;
            dList = [];
            alphaStacked = [];
            lambda = 0.5;
            initDepth = [];
            visibilityMatrix = [];
            angleList = [];
            termCriter = 10^-16;
            maxIter = 100;
            debugDisplay = true;
    end

    methods(Access = public)
            function [obj] = Depth(directionVectors, nng, depth, geodesic, nNeighbors, visibilityMatrix, angles)
                obj.directionVectors = directionVectors;
                obj.nng = nng;
                obj.depth = depth;
                obj.geodesic = geodesic;
                obj.nNeighbors = nNeighbors;
                obj.isDataSet = true;
                obj.MIN_DEPTH = obj.MAX_DEPTH/1000;
                obj.visibilityMatrix = visibilityMatrix;
                obj.angleList = angles;
            end

            function [obj] = setMaxDepth(obj, maxDepth)
                obj.MAX_DEPTH = maxDepth;
            end
            
            function [obj] = setInitDepth(obj, depth)
                obj.initDepth = depth;
            end
            
            function [obj] = setLambda(obj, lambda)
                obj.lambda = lambda;
            end
            
            function [obj] = updateGeodesic(obj, geodesic)
                obj.geodesic = geodesic;
            end
            
            function [obj] = setTerminationCriterion(obj, criter)
                obj.termCriter = criter;
            end
            
            function [obj] = setMaxIteration(obj, iter)
                obj.maxIter = iter;
            end
            
            function [obj] = switchOffDebugDisplay(obj)
                obj.debugDisplay = false;
            end

            function [depth, errorEvolution] = compute(obj)
                    assert(obj.isDataSet == true,'Depth@compute: set values of variables through constructor of <Depth> before calling this method.')

                    nPts = size(obj.directionVectors{1},2);            
                    nFiles =  size(obj.directionVectors,2);
                    nNieghbors = size(obj.nng, 2);
                    
                    prevError = realmax; currError = 0; innerIter = 0;
                    errorEvolution = [];
                        
                        while((innerIter < obj.maxIter ) && (abs(currError - prevError) > obj.termCriter))
                            for iFile = 1:nFiles
                                for iPt = 1:nPts         
                                    if(obj.visibilityMatrix(iFile, iPt) == 1)
                                        G = 0; H = 0; I = 0; N = 0; C = 0;
                                        nValidneighbors = 0;
                                        for iNeighbor = 1:nNieghbors
                                            indexQ = obj.nng(iPt, iNeighbor);
                                            
                                            if(obj.visibilityMatrix(iFile, indexQ) == 1)
                                                delQ =obj.depth{iFile}(indexQ);     
                                                d_jq = obj.geodesic{iPt, indexQ};

                                                prod = dot(obj.directionVectors{iFile}(:,iPt), obj.directionVectors{iFile}(:,indexQ));
%                                                 theta = obj.angleList{iFile}(iPt, iNeighbor); 
%                                                 prod = cos(theta); % previous option                                                

                                                D = (12 * delQ * prod);
                                                E = (delQ^2) - (d_jq)^2;

                                                G = G - D;          % H = H + ( (2 * E) + D^2 );  % ===> last version of H
                                                H = H + 2*( 2*delQ^2*(1 + (2*(prod^2)) ) - 2*d_jq^2);  
                                                I = I + ( (4*d_jq^2*delQ*prod) - (4*delQ^3*prod) );     N = N + E;
                                                C = C + delQ;
                                                
                                                nValidneighbors = nValidneighbors + 1;
                                            end
                                        end

                                        I = I - obj.lambda;
                                        F = 4 * nValidneighbors ;
                                        J = (1/4) * F;
                                        K = (1/3) * G;  L = (1/2) * H;
                                        M = I;  N = (2* N) - (obj.lambda * C);

                                        %[sol] = solution_general_cubic(F, G, H, I);
                                        [sol] = roots([F, G, H, I]);

                                        depthUpdate = -1;
                                        evalCost = realmax;                
                                        
                                        hasOneSol = false;

                                        for iSol = 1:numel(sol)
                                            if( (abs(imag(sol(iSol))) <10^(-30) ) && (sol(iSol) > 0))
                                                hasOneSol = true;
                                                [C] = obj.evaluateQuartic(J, K, L, M, N, real(sol(iSol)));
                                                [dD] = obj.evaluateDoubleDerivative(F, G, H, real(sol(iSol)));
                                                if( (C < evalCost) && (dD > 0) )
                                                    depthUpdate = real(sol(iSol));
                                                    evalCost = C;
                                                end
                                            end
                                        end
                                        
                                        if(~hasOneSol)
                                            warning('Did not find any relevant solution');
                                            fprintf('=====> Negative roots!\n');
                                        end

                                        if(depthUpdate > 0)
                                            obj.depth{iFile}(iPt) = depthUpdate;
                                        else
%                                             fprintf('**********************************************************************\n');
%                                             fprintf('**********************************************************************\n');
%                                             warning('Depth@compute: no real & positive roots found for the cubic! Details below:');
%                                             fprintf('File no.: %d,   feat no.: %d, inner AO iteration no.: %d\n', iFile, iPt, innerIter);
%                                             fprintf('Last known depth of this point = %d\n', obj.depth{iFile}(iPt));
%                                             fprintf('Cubic (derivated) polynomial coefficients: [%d, %d, %d, %d]\n', F, G, H, I);
%                                             fprintf('Quartic (main) polynomial coefficients: [%d, %d, %d, %d, %d]\n', J, K, L, M, N);
%                                             fprintf('Obtained solutions of cubic:\n');
%                                             for iP = 1:numel(sol)
%                                                 fprintf(' --------: Solution %d = %d\n', iP, sol(iP));
%                                             end
%                                             fprintf('**********************************************************************\n');
%                                             fprintf('**********************************************************************\n');
                                        end
                                    
                                    end
                                end
                            end
                            prevError = currError;
                            currError = evaluateTotalCostCombined(obj.directionVectors, obj.nng, obj.depth, obj.geodesic, obj.lambda, obj.visibilityMatrix);
                            errIdentifier = [currError 1];
                            if(innerIter == 0)
                                errorEvolution = errIdentifier;
                            else
                                errorEvolution = [errorEvolution; errIdentifier];
                            end
                            innerIter = innerIter + 1;
                            
                            if(obj.debugDisplay)
                                fprintf('iteration = %d, error = %d   [INNER]\n', innerIter, currError);
                            end
                        end
                        
                        depth = obj.depth;

            end
    end

    methods(Access = private)
        function [val] = evaluateQuartic(obj, a, b, c, d, e, X)
            val = (a .* X.^4) + (b .* X.^3) + (c .* X.^2) + (d .* X) + e;
        end       
        
        function [val] = evaluateDoubleDerivative(obj, F, G, H, delJ)
            val = (3 * F * delJ.^2) + (2 * G * delJ) + H;
        end
    end

end
