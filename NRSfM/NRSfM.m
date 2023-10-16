classdef NRSfM < handle
    %% Method for isometric NRSfM with Alternating Optimisation (AO) using 
    % Distance Constraint Function (DCF). The DCFs are solved by forming 
    % polynomials out of point pairs and the global minima is obtained by solving for 
    % the optimal geodesic distance using approximation.
    
    properties(Access = private)
        K = []; % [3x3] intrinsics matrix
        p = []; % keypoint correspondences
        pGTh = []; % GT keypoints
        v = []; % visibility matrix
        directionVectors = []; % derived direction vectors (unit normals) corresponding to each feature vector from correspondences
        Data = []; % the combined 'Data' struct loaded from file
        
        nK = 0; % no. of neighbors in the NNG
        terminationCriterion = 10^(-25); % termination criterion based on diff. of abs. error inbetween AO iterations
        innerTerminationCriterion = 10^(-25); % termination criterion of inner AO (depth) based on diff. of abs. error inbetween AO iterations
        maxIteration = 30; % max. permissible iteration for AO
        initialDepth = 1.0; % depth at which the 3D points are initialised at        
        maxInnerIteration = 30; % maximum no. of iteration permitted in the inner AO loop (for depth computation)        
        isDataAssigned = false; % flag indicating if data has been assigned, defaults to false
        isIntrinsicsAssigned = false; % flag indicating if intrinsics (K) has been assigned, defaults to false        
        bootstrapDepth = []; % to be used to update depth from other NRSfM method
        isBootstrapDepthAvailable = false;
        
        debugDisplay = true;
        
        geodesicsInit = [];
        isGeodesicInitializationProvided = false;
        
        lambda = 0;
        
        NNG = [];
        geodesics = [];
        depth = [];
        
        doneNRSfM = false
        
    end
    
    methods(Access = public)
      %% Constructor for the NRSfM class, sets the data 'struct'
        function [obj] = NRSfM(Data)
        % Constructor for the NRSfM class, sets the data 'struct'
        %
        % *PARAMETER*:
        %
        % _Data_: must be a 'struct' type variable with the following
        %               fields:
        %
        % * _Data_.p_: [1 x N] struct of [2 x M] arrays each
        %               containing feature correspondences
        % * _Data.Pgth_: [1xN] struct of [3 x M] arrays each
        %               containing GT as 3D points
        % * _Data.v_: [N x M] array of binary values (0 and 1
        %               only), indicating visibility status of each
        %               correspondences
        %
        % Where: N = no. of image files
        %               M = max. no. of correspondences in the
        %               sequences
        %
        % *Example:*
        %       To access the 3rd correspondence point
        %       of the 5th image, use:
        %       Data.p(5).p(:,3)   { for GT: Data.Pgth(5).P(:,3) }
        %       and this point is visible only if:  Data.v(5,3) == 1
        %            
            assert(isstruct(Data),'DataError: input data must be of type <struct>!');
            assert(size(Data.Pgth,2) > 1,'DataError: input GT size must be greater than one!');
            assert(size(Data.p,2) > 1,'DataError: input keypoint correspondence size must be greater than one!');
            assert(size(Data.v,2) > 1,'DataError: input visibility matrix size must be greater than one!');
            assert(size(Data.p,2) == size(Data.v,1),'DataError: number of files in correspondences must be equal to number of rows in visibility matrix');
            assert(size(Data.p,2) == size(Data.Pgth,2),'DataError: number of files in correspondences must be equal to number of files in GT');
                        
            sumV = sum(Data.v(1,:));
            assert(sumV == size(Data.v, 2), 'DataError: all keypoints in first image needs to be visible');
            
            obj.p = Data.p;
            obj.pGTh = Data.Pgth;
            obj.v = Data.v;
            obj.Data = cleanInvisiblePoints(Data);
            obj.isDataAssigned = true;
        end
        
      %% Set intrinsics for the NRSfM
        function [obj] = setIntrinsics(obj, intrinsics)
            % setIntrinsics: set intrinsics for the NRSfM
            %
            % *PARAMETER*:
            %    _K_ : [3 x 3] camera intrinsics matrix, where:
            %           K = [fx 0 cx; 0 fy cy; 0 0 1]
            %
            % *NOTE*: distortion parameters are not considered here, please undistort...
            %   ...images (and/or correspondences) before invoking this method
            %
            assert(size(intrinsics,1) == 3,'IntrinsicsError: K should be a [3x3] matrix');
            assert(size(intrinsics,2) == 3,'IntrinsicsError: K should be a [3x3] matrix');
            assert( ((intrinsics(1,1) > 0) && (intrinsics(2,2) > 0) && (intrinsics(1,3) > 0) && (intrinsics(2,3) > 0) ),'IntrinsicsError: K should be formatted as [fx 0 cx; 0 fy cy; 0 0 1]');
            
            obj.K = intrinsics;
            obj.isIntrinsicsAssigned = true;
        end
        %% Set termination criterion for AO
        function [obj] = setTerminationCriterion(obj, threshold)
            % setTerminationCriterion: set termination criterion for AO
            %
            obj.terminationCriterion = threshold;
        end
        
        %% Set termination criterion for inner AO loop for depth computation
        function [obj] = setInnerTerminationCriterion(obj, threshold)
            % setTerminationCriterion: set termination criterion for AO
            %
            obj.innerTerminationCriterion = threshold;
        end
        
        function [obj] = setInitialDepth(obj, depth)
            obj.bootstrapDepth = depth;
            obj.isBootstrapDepthAvailable = true;
        end
        
        function [obj] = setMaxIteration(obj, maxIter)
            obj.maxIteration = maxIter;
        end
        
        function [obj] = setInitialGeodesic(obj, initGeodesic)
            obj.geodesicsInit = initGeodesic;
            obj.isGeodesicInitializationProvided = true;
        end
        
        function [obj] = setLambda(obj, lambda)
            obj.lambda = lambda;
        end
        
        function [obj] = setMaxInnerIteration(obj, maxIter)
            obj.maxInnerIteration = maxIter;
        end
        
        function [obj] = switchOffDebugDisplay(obj)
            obj.debugDisplay = false;
        end
        
        %% The function to trigger NRSfM computation
        function [errorEvolution] = computeNRSfM(obj, numNeighbors)
            % computeNRSfM: the function to trigger NRSfM computation
            %
            if( obj.isDataAssigned && obj.isIntrinsicsAssigned )
                obj.nK = numNeighbors;
                [errorEvolution] = obj.doNRSfM();
            else
                error('MISSING INPUT: assign data and intrinsics first, before trying to invoke NRSfM\n');
            end
        end
        
        %% Recover reconstructed 3D points
        function [reconstruction] = get3DReconstruction(obj)
            for fIndex = 1:size(obj.v,1)
                Pest = obj.depth{fIndex}.*obj.directionVectors{fIndex};
                Pgt = obj.pGTh(fIndex).P;
                [scaledPoints] = scale3DPoints(Pest.', Pgt.');
                reconstruction{fIndex} = scaledPoints;
            end
        end
    end
    
    methods(Access = private)
       %% Actual implementation of the NRSfM
        function [errorEvolution] = doNRSfM(obj)
            
            % Setup
            [obj.directionVectors, depthSet] = getDirectionVectors(obj.Data, obj.K, obj.initialDepth, false);
            if(obj.isBootstrapDepthAvailable)
                depthSet = obj.bootstrapDepth;
            end
            [nng] = getNeighborhoodDuplicated(obj.Data, obj.nK);
            obj.NNG = nng;
            [angleList] = precomputeAngles(obj.directionVectors, obj.v, nng);
            
            if(obj.isGeodesicInitializationProvided)
                currError = evaluateTotalCostCombined(obj.directionVectors, nng, depthSet, obj.geodesicsInit, obj.lambda, obj.v);
            else
                currError = -1;
            end
            prevError = 0;   iteration = 0;   innerTermCriter = obj.terminationCriterion; errorEvolution = [];
            
            % Initialisation
            geodesicHandle = Geodesics(obj.directionVectors, nng, depthSet, obj.nK, obj.v, angleList);         
            depthHandle = Depth(obj.directionVectors, nng, depthSet, [], obj.nK, obj.v, angleList);               
            depthHandle.setLambda(obj.lambda);  depthHandle.setMaxIteration(obj.maxInnerIteration);
            depthHandle.setTerminationCriterion(obj.innerTerminationCriterion);
            
            if(obj.debugDisplay)
                timestmp = strrep(datestr(datetime('now')),' ','@');
                fprintf('\n\n<strong>Starting NRSfM at [%s]</strong>\n',timestmp);
            else
                [depthHandle] = depthHandle.switchOffDebugDisplay();
            end
            
            while( (abs(currError - prevError) > obj.terminationCriterion) && (iteration < obj.maxIteration))
                % Computing geodesic
                [geodesicHandle] = geodesicHandle.setDepth(depthSet);
                [geodesicSet] = geodesicHandle.compute();  
                
                prevError = currError;
                currError = evaluateTotalCostCombined(obj.directionVectors, nng, depthSet, geodesicSet, obj.lambda, obj.v);
                
                errIndex = [currError, 0];
                
                iteration = iteration + 1;  innerTermCriter = innerTermCriter/10;
                if(obj.debugDisplay)
                    fprintf(' ========> <strong>Iteration = %d</strong>, error = %d   [OUTER]\n', iteration, currError);                
                end

                % Computing depth
                depthHandle.updateGeodesic(geodesicSet);     
                depthHandle.setTerminationCriterion(innerTermCriter);
                [depthSet, innerErrorEvolution] = depthHandle.compute(); 
                
                errorEvolution = [errorEvolution; errIndex; innerErrorEvolution];

            end
            
            if(obj.debugDisplay)
                timestmp = strrep(datestr(datetime('now')),' ','_');
                fprintf('<strong>Ending NRSfM at [%s]</strong>\n\n',timestmp);
            end
            
            obj.geodesics = geodesicSet;
            obj.depth = depthSet;
            obj.doneNRSfM = true;
        end
    end
end