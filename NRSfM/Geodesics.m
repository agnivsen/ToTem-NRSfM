classdef Geodesics < handle
    
    properties(Access = private)
        directionVectors = [];
        nng = [];
        depth = [];
        geodesic = [] ;
        nNeighbors = 0;
        isDataSet = false;
        currAlpha = 0;
        currBeta = 0;
        shouldDebug = false;
        hasGT = false
        geodesicGT = [];
        MAX_DEPTH = 2.0;
        MIN_DEPTH = 0;
        initDepth = [];
        visibilityMatrix = [];
        angleList = [];
    end
    
    methods(Access = public)
        function [obj] = Geodesics(directionVectors, nng, depth, nNeighbors, visibilityMatrix, angles)
            obj.directionVectors = directionVectors;
            obj.nng = nng;
            obj.depth = depth;
            obj.nNeighbors = nNeighbors;
            obj.isDataSet = true;
            obj.MIN_DEPTH = obj.MAX_DEPTH/1000;
            obj.visibilityMatrix = visibilityMatrix;
            obj.angleList = angles;
        end
        
        function [obj] = debugModeOn(obj)
            obj.shouldDebug = true;
        end
        
        function [obj] = setDepth(obj, depth)
            obj.initDepth = depth;
        end
        
        function [obj] = setMaxDepth(obj, maxDepth)
            obj.MAX_DEPTH = maxDepth;
        end
        
        function [geodesic] = compute(obj)
            
            assert(obj.isDataSet == true,'Geodesics@compute: set values of variables through constructor of <Geodesics> before calling this method.')
            
            nPts = size(obj.directionVectors{1},2);
            
            for iAlpha = 1:nPts
                for iBeta = 1:obj.nNeighbors
                    if(obj.nng(iAlpha, iBeta) > 0)
                        
                        iNeigh = obj.nng(iAlpha, iBeta);
                        
                        obj.currAlpha = iAlpha;
                        obj.currBeta = iNeigh;
                        
                        nFiles =  size(obj.directionVectors,2);
                        
                        B = 0;
                        validFiles = 0;
                        
                        for iFile = 1:nFiles
                            if( ( obj.visibilityMatrix(iFile, iAlpha) == 1 ) && ( obj.visibilityMatrix(iFile, iNeigh) == 1 ) )
                                alpha = obj.depth{iFile}(1, iAlpha);
                                beta = obj.depth{iFile}(1, iNeigh);

%                                 prod = dot(obj.directionVectors{iFile}(:,iAlpha), obj.directionVectors{iFile}(:,iNeigh));
                                theta = obj.angleList{iFile}(iAlpha, iBeta);
                                prod = cos(theta); % previous option

                                B = B + ((alpha^2) + (beta^2) - (2 * alpha * beta * prod));
                                validFiles = validFiles + 1;
                            end
                        end
                        
                        obj.geodesic{iAlpha, iNeigh} = sqrt((1/validFiles) * B);
                        
                    else
                        obj.geodesic{iAlpha, iNeigh} = [];
                    end
                end
            end
            
            % % DO NOT normalize, this is unprincipled and makes the
            % % results worse!
            % [obj.geodesic] = obj.normalizeGeodesics(obj.geodesic);
            
            geodesic = obj.geodesic;
        end
    end
    
    methods(Access = private)
        function [geodesics] = normalizeGeodesics(obj, geodesicsInput)
            m = size(geodesicsInput,1);
            summedSquaredGeodesic = 0;
            geodesics = geodesicsInput;
            for ii = 1:m
                for jj = ii:m
                    if ~isempty(geodesicsInput{ii,jj})
                        summedSquaredGeodesic = summedSquaredGeodesic + geodesicsInput{ii,jj}^2;
                    end
                end
            end
            
            for ii = 1:m
                for jj = 1:m
                    geodesics{ii,jj} = geodesicsInput{ii,jj}/sqrt(summedSquaredGeodesic);
                end
            end
            
            
        end
    end
end