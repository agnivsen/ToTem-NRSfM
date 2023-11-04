% Copyright 2021 Fang Bai <fang dot bai at yahoo dot com>
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.


classdef DefGPA < handle
    
    
    % input data
    properties (Access = public)
        
        dim = 0
        
        PointClouds = struct ( 'Data', {}, 'Vis', {}, 'TransPrior', {} )
        
        numPointClouds = 0
        
        numPoints = 0
        
        dimFeatureSpace = 0
        
        warpModelFlag = 'AFFINE'  % TPS
        
        smoothParam = 0 % param for regularizations of warps
        
        centerDataAhead = false
        
        flagJointlyOptimizeLamdaRigidTransformation = true
        
        verbose = true

    end
    
    
    % output properties
    properties (Access = public)
        
        mShape = []
        
        mWeights = {}
        
        mPoses = {}
        
        mSqrtLambda = []
        
        ReferenceCovariancePrior = []
        
        TPS = struct ( 'IntrisicMatr', {}, 'ControlPoints', {}, 'InternalSmoothParam', {}, 'Weights', {})
        
        AFFINE = struct ('Affine', {}, 'Translation', {});
        
        FeatureClouds = struct ( 'FeatureMatr', {}, 'Regularization', {} )
        
        P = []  % P-Matrix

        rsd_ref = 0
        
        rsd_dat = 0
        
        rigidity = 0
        
        time = []
        
    end
    
    
    % internal private properties
    properties (Access = private)
        
        maxNum = 1000000  % the celing number, should never be reached.
        
        EuclideanTransformationSign = []
        
        fileExt = 'txt'
              
    end
    
    
    
    methods (Access = public)
 
        % class constructor
        function this = DefGPA ( modelFlag )            
            this.warpModelFlag = 'AFFINE';            
            if (nargin > 0)            
                this.warpModelFlag = modelFlag;            
            end            
        end
        
               

        function addDataByArray(this, DataCell, VisVecCell)
            for ii = 1 : size(DataCell,2)
                this.addPointCloud (DataCell{ii}, VisVecCell{ii});
            end
            if this.dim == 2
                this.smoothParam = 10;
            end
            if this.dim == 3
                this.smoothParam = 0.01;
            end
        end
        
        
        function run (this, smoothParam, flagJointlyOptimizeLamdaRigidTransformation)            
            if exist('smoothParam', 'var') && ~isempty(smoothParam)
                this.smoothParam = smoothParam;
            end
            if exist('flagJointlyOptimizeLamdaRigidTransformation', 'var') && ~isempty(flagJointlyOptimizeLamdaRigidTransformation)
                this.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
            end 
             this.RunFullProcedure();
        end
        
        
        function arrayTrainingStatistics = OptimizeSmoothParam (this, params,  Nfold)
            if ~exist('params', 'var') || (exist('params', 'var') && isempty(params))
                if this.dim == 3
                    params = [0.001, 0.01, 0.1, 1];
                end
                if this.dim == 2
                    params = [0.1, 1, 10, 100];
                end
            end
            if ~exist('Nfold', 'var') || (exist('Nfold', 'var') && isempty(Nfold))
                Nfold = 5;
            end
            arrayTrainingStatistics.SmoothParam = zeros(1, length(params));
            arrayTrainingStatistics.CVE = zeros(1, length(params));
            arrayTrainingStatistics.RMSE_ref = zeros(1, length(params));
            arrayTrainingStatistics.RMSE_dat = zeros(1, length(params));
            for ii = 1 : length(params)
                smooth_param = params(ii);
                this.smoothParam = smooth_param;
                this.RunFullProcedure();
                [RMSE_ref, RMSE_dat] = this.rsdError();
                cve = this.run_N_Fold_CrossValidation (Nfold);
                arrayTrainingStatistics.SmoothParam(ii) = smooth_param;
                arrayTrainingStatistics.CVE(ii) = cve;
                arrayTrainingStatistics.RMSE_ref(ii) = RMSE_ref;
                arrayTrainingStatistics.RMSE_dat(ii) = RMSE_dat;
            end
        end


        
        function [cve, predictedPointClouds] = run_N_Fold_CrossValidation (this, N)
            this.RunFullProcedure();
            gaugeShape = this.mShape;
            predictedPointClouds = cell(1, this.numPointClouds);
            for ii = 1 : this.numPointClouds
                predictedPointClouds{ii} = 0 * this.PointClouds(ii).Data;
            end
            N = min(this.numPoints, N);
            for ii = 1 : N
                vtmp = (ii : N : this.numPoints);
                testingPointsFlag = false(1, this.numPoints);
                testingPointsFlag(vtmp) = true;
                trainingPointsFlag = ~logical(testingPointsFlag);
                % compute the transformed test points
                [transform_points_cell, transform_vis_cell, refShape] = runTrainingTestData(this, trainingPointsFlag, testingPointsFlag);
                % correct gauge and fill the predicted points
                [sR, t] = this.SimilarityProcrustes (refShape, gaugeShape(:, trainingPointsFlag));
                for kk = 1 : this.numPointClouds
                    predictedPointClouds{kk}(:, testingPointsFlag) = (sR * transform_points_cell{kk} + t) .* transform_vis_cell{kk};
                end
            end

            % compute the mean template
            if(0)                
                ref_point_cloud = zeros(this.dim, this.numPoints);
                ref_num_vis = zeros(1, this.numPoints);
                for kk = 1 : this.numPointClouds
                    pts = predictedPointClouds{kk} .* this.PointClouds(kk).Vis;
                    ref_point_cloud = ref_point_cloud + pts;
                    ref_num_vis = ref_num_vis + this.PointClouds(kk).Vis;
                end
                for jj = 1 : this.numPoints
                    if (ref_num_vis(jj) > 0)
                        ref_point_cloud(:, jj) = ref_point_cloud(:, jj) ./ ref_num_vis(jj);
                    end
                end
            else
                ref_point_cloud = gaugeShape;
            end
            % CVE is computed using the discrepancy between the predicted-point-clouds and the global reference shape
            normalizer = 0; cve = 0;
            for ii = 1 : this.numPointClouds
                vis =  this.PointClouds(ii).Vis;
                pointsDeviation = (predictedPointClouds{ii} - ref_point_cloud) .* vis;
                sqrtError =  norm(pointsDeviation, 'fro');
                cve = cve + sqrtError * sqrtError;
                normalizer = normalizer + nnz(vis);                
            end
            cve = sqrt(cve/normalizer);
            % reset the estimate to the status of using all the points.
            this.RunFullProcedure(); 
%             for ii = 1 : this.numPointClouds
%                 VisVecCell{ii} = this.PointClouds(ii).Vis;
%             end
%             cve = this.PointCloudsVariationError (predictedPointClouds, VisVecCell);
        end


        
        function addPointCloud (this, dataMatr, visVec )
            cnt = this.numPointClouds + 1;
            this.PointClouds(cnt).Data = double(dataMatr);            
            this.PointClouds(cnt).Vis = ( sum (abs(dataMatr) < this.maxNum) == size(dataMatr, 1) );            
            if (nargin == 3)                
                this.PointClouds(cnt).Vis = logical(visVec);                
            end            
            this.numPointClouds = cnt;            
            this.dim = size(dataMatr, 1);            
            this.numPoints = size(dataMatr, 2);            
        end


        function pImageArray = transformPoints (this, pValArray, cloudIndex )            
            pImageArray = pValArray + this.PointClouds(cloudIndex).TransPrior;             
            for jj = 1 : size(pImageArray, 2)                
                pImageArray(:,jj) = this.mWeights{cloudIndex}' *  this.fMapping (pImageArray(:,jj), cloudIndex );                
            end            
        end
        
        function pImageArray = affineTransformPoints (this, pValArray, cloudIndex )            
            pImageArray = pValArray + this.PointClouds(cloudIndex).TransPrior;             
        end
        
        
        function ptCloud = generatePointCloud (this, cloudIndex)            
            ptCloud = this.inverseWarpMapping (cloudIndex, this.mShape) - this.PointClouds(cloudIndex).TransPrior;            
        end
        
        
        function [rsd_ref, rsd_dat] = rsdError (this)            
            rsd_ref = 0; rsd_dat = 0; v = 0;            
            for ii = 1 : this.numPointClouds                
                v = v + nnz(this.PointClouds(ii).Vis);                
                % reference space cost                
                Diff = this.mShape - this.mWeights{ii}' * this.FeatureClouds(ii).FeatureMatr;                
                dr = norm( Diff .* this.PointClouds(ii).Vis, 'fro');
                rsd_ref = rsd_ref + dr * dr;                
                % datum space cost                
                Diff = this.inverseWarpMapping (ii, this.mShape) - this.PointClouds(ii).Data;
                dd = norm( Diff .* this.PointClouds(ii).Vis, 'fro');                
                rsd_dat = rsd_dat + dd * dd;                
            end                        
            rsd_ref = sqrt(rsd_ref/v);            
            rsd_dat = sqrt(rsd_dat/v);            
            this.rsd_ref = rsd_ref;  this.rsd_dat = rsd_dat;
        end
        
        
        function rig = rigidityScore (this) 
            rig = zeros(this.numPointClouds, 1);            
            for ii = 1 : this.numPointClouds                
                D = this.PointClouds(ii).Data .* this.PointClouds(ii).Vis;                
                Dtrans = ( this.mWeights{ii}' * this.FeatureClouds(ii).FeatureMatr )  .* this.PointClouds(ii).Vis;                 
                rig (ii) = this.DistScore (D(:, logical(this.PointClouds(ii).Vis)), Dtrans(:, logical(this.PointClouds(ii).Vis)));                
            end            
            this.rigidity = 1 - rig;            
        end
        

        function readDir (this, dirPath)            
            fileHandles = dir ([dirPath, '/*.', this.fileExt]);           
            dataMatr = [];                        
            for ii = 1 : size (fileHandles, 1)                
                fid = fopen ([dirPath, '/', fileHandles(ii).name], 'r');                
                lstr = fgetl (fid);                
                pcnt = 0;                
                while ischar (lstr)                    
                    if (size (lstr, 2) > 0)                        
                        pcnt = pcnt + 1;                        
                        dataMatr(:, pcnt) = str2double(split(lstr));                        
                    end                    
                    lstr = fgetl (fid);                    
                end                
                this.addPointCloud (dataMatr);                
                fclose (fid);                
            end            
        end
        
        
        function writeShape (this, outputFile)            
            fid = fopen ([outputFile, '.xyz'], 'w');
            formatSpec = '%f %f %f\n';            
            for ii = 1 : this.numPoints                
                v = this.mShape(:, ii);                
                fprintf (fid, formatSpec, v);           
            end            
            fclose (fid);            
        end

    end
    
    
    
    

    
    methods (Access = private)
        

        function RunFullProcedure(this)
            tic;

            if ~this.flagJointlyOptimizeLamdaRigidTransformation
                this.calReferenceCovariancePrior('ascend'); % ascending order
            else
                this.ReferenceCovariancePrior = ones(this.dim, 1);
                transPrior = zeros(this.dim, 1);
                for ii = 1 : this.numPointClouds
                    this.PointClouds(ii).TransPrior = transPrior * logical(this.centerDataAhead);
                end
            end


            this.initWarpModel ();

            this.toFeatureCloud();

            this.setP(); % P 1 = 0

            if norm(this.P * ones(length(this.P), 1)) > 1e-6
                fprintf(2, 'Fatal Error: P*1 = 0 is violated. ||P*1|| = %.15f \n', norm(this.P * ones(length(this.P), 1), 'fro'));
                %return;
            end

            this.P = this.P + (this.numPointClouds/this.numPoints); % soft regularization

            [EV, eigval, ~] = this.byEigen (this.P, this.dim, 'smallest'); % descending order

            this.mShape = EV';

            if ~this.flagJointlyOptimizeLamdaRigidTransformation
                this.eliminateReflection();
            end

            this.setWeight();


            if this.flagJointlyOptimizeLamdaRigidTransformation
                this.InitializeLambda();
                rigid_cost = this.InitializePoses();
                if (this.verbose)
                    fprintf(1, 'iter[0]   cost = %f  sqrtLambda = (', rigid_cost);  fprintf(1, ' %f ', this.mSqrtLambda); fprintf(1, ')\n');
                end
                for ii = 1 : 100
                    [norm_update, rigid_cost] = this.PoseExtractionByNonlinearLeastSquaresOneIteration();
                    if (this.verbose)
                        fprintf(1, 'iter[%d]  cost = %f  ||update|| = %f  sqrtLambda = (', ii, rigid_cost, norm_update);  fprintf(1, ' %f ', this.mSqrtLambda); fprintf(1, ')\n');
                    end
                    if norm_update < 1e-7
                        break;
                    end
                end
                this.ReferenceCovariancePrior = this.mSqrtLambda .* this.mSqrtLambda;
                this.mShape = diag(this.mSqrtLambda) * this.mShape;
                this.setWeight();
            else
                this.mSqrtLambda = sqrt(this.ReferenceCovariancePrior);
                this.mShape = diag(this.mSqrtLambda) * this.mShape;
                this.setWeight();
                this.ConstructOptimalRigidTransform();
            end

            this.time(1) = toc;

        end

        
        function [transform_points_cell, transform_vis_cell, refShape] = runTrainingTestData(this, trainingPointsFlag, testingPointsFlag)
            AllPointClouds = this.PointClouds;
            AllNumPoints = this.numPoints;
            this.PointClouds = struct ( 'Data', {}, 'Vis', {}, 'TransPrior', {});
            this.numPointClouds = 0;
            for ii = 1 : length(AllPointClouds)
                data = AllPointClouds(ii).Data(:,  logical(trainingPointsFlag));
                vis = AllPointClouds(ii).Vis(:,  logical(trainingPointsFlag));
                this.addPointCloud (data, vis);
            end
            this.RunFullProcedure();
            % CVE: by testing points
            for ii = 1 : this.numPointClouds
                test_pts = AllPointClouds(ii).Data(:,  logical(testingPointsFlag));
                test_vis = AllPointClouds(ii).Vis(:,  logical(testingPointsFlag));
                transform_points_cell{ii} = this.transformPoints(test_pts, ii);
                transform_vis_cell{ii} = test_vis;
            end
            refShape = this.mShape;
            this.PointClouds = AllPointClouds;
            this.numPoints = AllNumPoints;
        end
        

        % It is better to calculate this cost using the estimated reference shape
        % The optimal R and t is extracted from the nonlinear transformation.
        % || R * Data + t - TransformShape || can keep reducing for big mu
        % || R * Data + t - RefShape || will reduce slower for big mu, as it's not minimized
        function cost = RigidCost (this)
            % compute the initial cost
            cost = 0;
            num_vis = 0;
            for ii = 1 : this.numPointClouds
                num_vis = num_vis + nnz(this.PointClouds(ii).Vis);
                R = this.mPoses{ii}(:, 1:this.dim);
                t = this.mPoses{ii}(:, this.dim+1);
                %Diff = (R * this.PointClouds(ii).Data + t) - this.mWeights{ii}' * this.FeatureClouds(ii).FeatureMatr;
                Diff = (R * this.PointClouds(ii).Data + t) - this.mShape;
                normtmp = norm( Diff .* this.PointClouds(ii).Vis, 'fro');
                cost = cost + normtmp * normtmp;
            end
            cost = sqrt(cost/num_vis);
        end
        
        
        function ConstructOptimalRigidTransform (this)
            for ii = 1 : this.numPointClouds
                Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));
                TransPt = full(this.mWeights{ii}' * this.FeatureClouds(ii).FeatureMatr(:, logical(this.PointClouds(ii).Vis)));
                meanTransPt = mean(TransPt, 2);  meanPt = mean(Pt, 2);
                R = (Pt - meanPt) * (TransPt - meanTransPt)';
                [U, ~, V] = svd(R); det_VU = det(V*U'); % this has the same sign as det(R). If det(R) > 0, then det(V'U) > 0
                diagS = eye(this.dim);  diagS(end, end) = det_VU;
                R = V * diagS * U';
                if (det_VU < 0)
                    fprintf(2, 'Warning @Det(VU) < 0 in orthogonal Procrustes \n');
                end
                t = - (R * meanPt - meanTransPt);
                this.mPoses{ii} = [R, t];
            end
        end
        
    end
    
    

    

    % functions used for Monte-Carlo sampling
    methods (Access = private, Static = true)
        
        function [idx_opt, idx_left_right] = ProcessMetricCost (metric_cost)
            [~, idx_opt] = min(metric_cost);
            idx_left_right = [idx_opt-1,  idx_opt+1];
            if (idx_opt == 1)
                idx_left_right(1) = 1;
            end
            if (idx_opt == length(metric_cost))
                idx_left_right(2) = length(metric_cost);
            end
        end
        
        
        function sample_array = SampleRangeEvenely (range_val1,  range_val2, num_samples)
            if (range_val1 < range_val2)
                sample_array = range_val1 :  ((range_val2 - range_val1)/(num_samples-1)) : range_val2;
            else
                sample_array = range_val2 :  ((range_val1 - range_val2)/(num_samples-1)) : range_val1;
            end
        end
        
    end
    
    
    
    % methods used to constructed internal matrices
    methods (Access = private)
  
        function toFeatureCloud (this)
            
            if strcmp (this.warpModelFlag, 'AFFINE')
                
                this.dimFeatureSpace = this.dim + 1;
                
            end
            
            for ii = 1 : this.numPointClouds
                
                this.FeatureClouds(ii).FeatureMatr = zeros ( this.dimFeatureSpace, this.numPoints );
                
                for jj = 1 : this.numPoints
                    
                    this.FeatureClouds(ii).FeatureMatr(:, jj) = this.fMapping ( this.PointClouds(ii).Data(:, jj), ii );
                    
                end
                
                this.FeatureClouds(ii).FeatureMatr = this.FeatureClouds(ii).FeatureMatr .* this.PointClouds(ii).Vis;
                
                this.FeatureClouds(ii).Regularization = this.getRegularization (ii);
                
            end
            
        end
        
        
        
        
        function setP (this)
            
            Q = zeros(this.numPoints, this.numPoints);
                        
            for ii = 1 : this.numPointClouds
                
                BMatr = this.FeatureClouds(ii).FeatureMatr;
                
                ZZMatr = this.FeatureClouds(ii).Regularization;

                Q = Q + BMatr' * pinv (full(BMatr * BMatr' + ZZMatr)) * BMatr;
                
            end
            
            vis = zeros(1, this.numPoints);
            
            for ii = 1 : this.numPointClouds
                
                vis = vis + this.PointClouds(ii).Vis;
                
            end
            
            this.P = diag(vis) - Q;
            
        end
        
        
        function v = calReferenceCovariancePrior (this, specificStr)
            
            orderingDirection = 'ascend';
            
            if nargin == 2
                
                orderingDirection = specificStr;
                
            end
            
            sumVis = zeros(1, this.numPoints);
            
            for ii = 1 : this.numPointClouds
                 
                sumVis = sumVis + this.PointClouds(ii).Vis;
            
            end
            
            sumVis = 1 ./ sumVis;
            
            eigVal = zeros (this.dim, this.numPointClouds);
            
            for ii = 1 : this.numPointClouds 
                 
                % Partial Shape: complete a full shape by point occurrence in other shapes
                if ( sum(this.PointClouds(ii).Vis) ~= this.numPoints ) 
                    
                    sumD = zeros(this.dim, this.numPoints);
                    
                    for jj = 1 : this.numPointClouds
                        
                        visFlag = logical(this.PointClouds(ii).Vis .* this.PointClouds(jj).Vis);
                        
                        [sRij, tij] = this.SimilarityProcrustes ( full(this.PointClouds(jj).Data(:, visFlag)),  full(this.PointClouds(ii).Data(:, visFlag)) );
                        
                        sumD = sumD + (sRij * this.PointClouds(jj).Data.*this.PointClouds(jj).Vis  + tij);
                         
                        if det(sRij) < 0
                            fprintf(2, '\n det(SRij ) = %f, by i = %d, j = %d\n', det(sRij), ii, jj); 
                        end
                        
                    end
                    
                    fullD = this.PointClouds(ii).Data .* this.PointClouds(ii).Vis + sumD .*  ( sumVis .* (1 - this.PointClouds(ii).Vis) );
                    
                else % Full Shape
                    
                    fullD = this.PointClouds(ii).Data;
                    
                end
                
                transPrior = - mean(fullD, 2);
                
                fullD = fullD + transPrior;
                
                tmp = fullD * fullD';
                
                eigVal(:, ii) = sort ( eig ( (tmp+tmp')/2 ), orderingDirection );
                
                % the following code applies a prior translation to all data point-clouds
                % it's not necessary to center the data ahead if the transformation model can handle translation.
                % Otherwise center the data ahead to compensate the translation.
                this.PointClouds(ii).TransPrior = transPrior * logical(this.centerDataAhead); 
                
            end
            
            
            for ii = 1 : this.numPointClouds
                
                this.PointClouds(ii).Data = this.PointClouds(ii).Data + this.PointClouds(ii).TransPrior;
                
                this.PointClouds(ii).Data = this.PointClouds(ii).Data .* this.PointClouds(ii).Vis;
                
            end
            
            
            scaleFactor = sqrt(sum(eigVal, 1));
            
            UnitSingularVal = sqrt (eigVal) ./ scaleFactor;
            
            [U, ~, ~] = svd( UnitSingularVal, 'econ');
            
            v = U(:, 1) * mean(scaleFactor);
            
            this.ReferenceCovariancePrior = v .* v;
            
        end
        
        
        
        function eliminateReflection (this)
            
            vis = logical (this.PointClouds(1).Vis);
            
            D1 = this.PointClouds(1).Data(:, vis);
                       
            R = this.OrhoProcurestes (D1, this.mShape(:, vis)); % will be zero-centered inside the function.
            
            if (det(R) < 0)
            
                this.mShape(1,:) = - this.mShape(1,:); % flip the sign of the first row.
                
                R = this.OrhoProcurestes (D1, this.mShape(:, vis));
                
                if (det(R) < 0)
                    fprintf(2, '\n Error Occurred in Eliminating Reflection! \n');
                end
                             
            end
            
        end
        
        
        
    end
    
    
  
    

    
    
    % utility functions
    methods (Access = private, Static = true)
        
        % a function used to construct a sparse matrix by given sublocks
        % this function put a matrix A at the position (row, col) of a sparse matrix.
        % use " sparse (i, j, v) " to check the resulting sparse matrix
        function [i, j, v] = putMatrixBlock (A, row, col)
            
            [m, n] = size(A);
            
            v = full (reshape (A, 1, m * n));
            
            i = kron (ones(1,n), row:row+m-1);
            
            j = kron (col:col+n-1, ones(1,m));
            
        end
        

        % distance score between two shapes
        % used to evaluate rigidity of transformations
        function distScore = DistScore (D, Dtrans)
            
            AS = alphaShape (D');
            TR = triangulation(AS.alphaTriangulation(), AS.Points);
            edges = TR.edges();
            
            distD = 0; distDiff = 0;
            
            for jj = 1 : length(edges)
                
                e = edges(jj, :);
                
                dist = norm( D(:, e(1)) - D(:, e(2)) );
                
                distErr = abs( norm( D(:, e(1)) - D(:, e(2)) ) - norm ( Dtrans(:, e(1)) - Dtrans(:, e(2)) ) );
                
                distD = distD + dist; distDiff = distDiff + distErr;                
                
            end
            
            distScore = distDiff/distD;
                        
        end
        
    end
    
    
    
    % static methods of solvers
    methods (Access = private, Static = true)
        
        % K, the matrix to decompose,
        % d, the number of eigen vectors required.
        function [V, eigValEigen, eigenValueAll] = byEigen (K, d, str)
            
            offset = 0;

            [Q, D ] = eig( (K+K')/2 );
            
            eigVal = diag(D);
            
            [eigVal, Index] = sort(eigVal, 'descend');  Q = Q (:, Index);

            if  strcmp (str, 'largest')
                
                V = Q (:, 1+offset:d+offset);
                
                eigValEigen =  eigVal (1+offset:d+offset);

            end
            
            if strcmp (str, 'smallest')
                
                V = Q (:, end-d+1-offset:end-offset);
                
                eigValEigen =  eigVal (end-d+1-offset:end-offset);

            end
            
            eigenValueAll = eigVal;
            
%             if  norm( V' * V - eye(d), 'fro' ) > 1e-12
%                 
%                 fprintf(2, 'Eigen vectors are not orthnormal!.');
%                 
%             end

        end
        
        

        % The solution to the similarity Procrustes problem
        % sR, t = minimize || sR * D1  + t  - D2 ||
        function [sR, t] = SimilarityProcrustes (D1, D2)
            
            meanD1 = mean(D1, 2);
            
            meanD2 = mean(D2, 2);
            
            M = (D1 - meanD1) * (D2 - meanD2)';
            
            N = (D1 - meanD1) * (D1 - meanD1)';
            
            [U, ~, V] = svd(M);
            
            R = V * U';            
            
            s = trace(R * M) / trace(N);
            
            sR = s * R;
            
            t =  meanD2 - sR * meanD1;
            
        end
        
        

        % The solution to the Euclidean Procrustes problem
        % R, t = minimize || R * D1  + t  - D2 ||
        function [R, t] = EuclideanProcrustes (D1, D2)

            meanD1 = mean(D1, 2);

            meanD2 = mean(D2, 2);

            M = (D1 - meanD1) * (D2 - meanD2)';

            N = (D1 - meanD1) * (D1 - meanD1)';

            [U, ~, V] = svd(M);

            R = V * U';

            t =  meanD2 - R * meanD1;

        end

        

        % The solution to the Orthogonal Procrustes problem
        % R = minimize || R * D1  - D2 ||
        function R = OrhoProcurestes (D1, D2)
            
            meanD1 = mean(D1, 2);
            
            meanD2 = mean(D2, 2);
            
            M = (D1 - meanD1) * (D2 - meanD2)';
            
            [U, ~, V] = svd(M);
            
            R = V * U';
            
        end
        
                
    end

    
    % Methods that has to be modified for new warp models
    % Methods to be overloaded based on different warp models.
    % implement with flags as an option
    methods (Access = private)
        
        function setWeight (this)
            
            for ii = 1 : this.numPointClouds
                
                BMatr = this.FeatureClouds(ii).FeatureMatr;
                
                ZZMatr = this.FeatureClouds(ii).Regularization;
                
                this.mWeights{ii} = pinv (full(BMatr * BMatr' + ZZMatr)) * full(BMatr * this.mShape');
                
                if strcmp (this.warpModelFlag, 'TPS')
                    
                    this.TPS(ii).Weights = this.mWeights{ii};
                    
                end
                
                if strcmp (this.warpModelFlag, 'AFFINE')
                                        
                    this.AFFINE(ii).Affine = this.mWeights{ii}(1:end-1, :)';
                    
                    this.AFFINE(ii).Translation = this.mWeights{ii}(end, :)';
                    
                end
                
                if strcmp (this.warpModelFlag, 'your new warp model')

                end
                
            end
            

        end
        
        
        % functions used to intialize constants to be used later on
        function initWarpModel (this)
            
            if strcmp (this.warpModelFlag, 'TPS')
                
                for ii = 1 : this.numPointClouds
                    
                    this.TPS(ii).ControlPoints =  ThinPlateSpline.setControlPoints (this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)), this.dimFeatureSpace);
                    
                    this.TPS(ii).InternalSmoothParam = 0; % diagonal elements of TPS kernel matrix. default = 0
                    
                    this.TPS(ii).IntrisicMatr = ThinPlateSpline.setTPSIntrisicMatr (this.TPS(ii).ControlPoints, this.TPS(ii).InternalSmoothParam);
                    
                end
                
                return
                
            end
            
            if strcmp (this.warpModelFlag, 'AFFINE')
                
                return;
                
            end
            
            if strcmp (this.warpModelFlag, 'your new warp model')
                
                return;
                
            end
            
        end
        
        
        
        % mapping R^3 -> R^l : point to feature space
        function fImage = fMapping (this, pVal, cloudIndex)
           
            if strcmp (this.warpModelFlag, 'TPS')
                
                nlv = ThinPlateSpline.getEtaTPS (this.TPS(cloudIndex).ControlPoints, pVal);

                fImage = this.TPS(cloudIndex).IntrisicMatr' * nlv;
                
                return
                
            end
            
            if strcmp (this.warpModelFlag, 'AFFINE')
                
                fImage = [pVal; 1];
                
                return;
                
            end
            
            if strcmp (this.warpModelFlag, 'your new warp model')
                
                return;
                
            end
            
        end
        
        
        % set the regularization for each FeatureCloud
        function M = getRegularization (this, cloudIndex)
            
            if strcmp (this.warpModelFlag, 'TPS')
                
                M = 8 * pi * this.TPS(cloudIndex).IntrisicMatr (1:this.dimFeatureSpace, :);
                
                M = (this.smoothParam * nnz(this.PointClouds(cloudIndex).Vis) ) * M;
                
                return
                
            end
            
            if strcmp (this.warpModelFlag, 'AFFINE')
                
                M = zeros(this.dimFeatureSpace);
                
                return;
                
            end
            
            if strcmp (this.warpModelFlag, 'your new warp model')
                
                return;
                
            end
            
        end
        
        
        % implement this methods if you want to compute the data space cost
        function inverseWarpCloud = inverseWarpMapping (this, cloudIndex, inShape)

            inverseWarpCloud = zeros (size(inShape, 1), size(inShape, 2));
            
            if strcmp (this.warpModelFlag, 'TPS')
                
                M = this.TPS(cloudIndex).ControlPoints * (ThinPlateSpline.setTPSIntrisicMatr(this.TPS(cloudIndex).Weights', this.TPS(cloudIndex).InternalSmoothParam))';
                
                for jj = 1 : this.numPoints
                    
                    inverseWarpCloud(:,jj) = M * ThinPlateSpline.getEtaTPS(this.TPS(cloudIndex).Weights', inShape(:,jj));
                    
                end
                
                return;
                
            end
            
            if strcmp (this.warpModelFlag, 'AFFINE')
                
                inverseWarpCloud = inv(this.AFFINE(cloudIndex).Affine) * ( inShape - this.AFFINE(cloudIndex).Translation);

                return;
                
            end
            
            if strcmp (this.warpModelFlag, 'your new warp model')
                
                return;
                
            end
                
        end
        
        
        
        
    end
    
    
    
    
    
    properties (Access = private)
        
        TmpA_ijv = struct('i_vec', [], 'j_vec', [], 'v_vec', []);
        
    end
    
      
    % methods used to compute Lambda and poses
    methods (Access = private)
        
        
        
        function cost = InitializePoses (this)
            if abs(sum(this.EuclideanTransformationSign)) ~= this.numPointClouds
                fprintf(2, 'Consistent optimal rigid transformation to the reference point-cloud DO NOT exist. \n');
                sign_of_euclidean_transformations = this.EuclideanTransformationSign
            end
            for ii = 1 : this.numPointClouds
                Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));
                St = full(this.mWeights{ii}' * this.FeatureClouds(ii).FeatureMatr(:, logical(this.PointClouds(ii).Vis)));
                meanSt = mean(St, 2);
                meanPt = mean(Pt, 2);
                R = (Pt - meanPt) * (St - meanSt)' * diag(this.mSqrtLambda);
                [U, ~, V] = svd(R); det_VU = det(V*U'); % this has the same sign as det(R). If det(R) > 0, then det(V'U) > 0
                diagS = eye(this.dim);  diagS(end, end) = det_VU;
                R = V * diagS * U';
%                 if (det_VU < 0)
%                     fprintf(2, 'Warning @Det(VU) < 0 in orthogonal Procrustes \n');
%                 end
                t = - (R * meanPt - diag(this.mSqrtLambda) * meanSt);
                this.mPoses{ii} = [R, t];
            end
            % compute the initial cost
            cost = 0;
            num_vis = 0;
            for ii = 1 : this.numPointClouds
                num_vis = num_vis + nnz(this.PointClouds(ii).Vis);
                Pt = this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis));
                St = this.mWeights{ii}' * this.FeatureClouds(ii).FeatureMatr(:, logical(this.PointClouds(ii).Vis));
                R = this.mPoses{ii}(:, 1:this.dim);
                t = this.mPoses{ii}(:, this.dim+1);
                vrhs = (R * Pt + t) - diag(this.mSqrtLambda) * St;
                normtmp = norm(vrhs, 'fro');
                cost = cost + normtmp * normtmp;
            end
            cost = sqrt(cost/num_vis);
        end
        


        function Lambda = InitializeLambda (this)
            this.EuclideanTransformationSign = zeros(1, this.numPointClouds);
            RgLambdaRg = zeros(this.dim ,this.dim);
            for ii = 1 : this.numPointClouds
                Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));
                Qt = this.mShape(:, logical(this.PointClouds(ii).Vis));
                %Qt = full(this.transformPoints(Pt, ii));
                zcQt = Qt - mean(Qt, 2);
                zcPt = Pt - mean(Pt, 2);
                invLt = full(zcPt * zcQt' * pinv(full(zcQt * zcQt')));
                RgLambdaRg = RgLambdaRg + invLt' * invLt;
                this.EuclideanTransformationSign(ii) = sign(det(invLt));  % the sign of affine transformation                
            end

            RgLambdaRg = RgLambdaRg ./ this.numPointClouds;
            [Rg, diagLambda] = eig( (RgLambdaRg + RgLambdaRg')/2 );          

            [Lambda, sindex] = sort(diag(diagLambda), 'ascend');
            Rg = Rg(:, sindex); Rg = Rg';

            %testouput_lambda = Lambda
            if (nnz(Lambda < 0))
                fprintf(2, 'Error@ There exist negative components in Lambda\n');
            end
            this.mSqrtLambda = sqrt(Lambda);

            if (sum(this.EuclideanTransformationSign < 0) > sum(this.EuclideanTransformationSign > 0))
                this.mShape(1,:) = - this.mShape(1,:); % flip the sign of the first row.
                this.EuclideanTransformationSign = - this.EuclideanTransformationSign;
            end
        end


        % reflections doesn't affect the estimate of Lambda
%         function Lambda = InitializeLambda (this)
%             this.EuclideanTransformationSign = zeros(1, this.numPointClouds);
%             AMatrix = zeros(this.dim * this.dim * this.numPointClouds, this.dim);
%             bVector = zeros(this.dim * this.dim * this.numPointClouds, 1);
%             n_elements = this.dim * this.dim;
%             for ii = 1 : this.numPointClouds
%                 Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));
%                 St = full(this.mWeights{ii}' * this.FeatureClouds(ii).FeatureMatr(:, logical(this.PointClouds(ii).Vis)));
%                 zcSt = St - mean(St, 2);
%                 zcPt = Pt - mean(Pt, 2);
%                 Lt = zcSt * zcPt' * pinv(zcPt * zcPt');
%                 this.EuclideanTransformationSign(ii) = sign(det(Lt));  % the sign of affine transformation
%                 %  || L' * Lambda * L - Identity ||_F = ||A * lambda - b||_2
%                 if (length(Lt) == 2)
%                     H1 = Lt(1,:)' * Lt(1,:);
%                     H2 = Lt(2,:)' * Lt(2,:);
%                     A = [reshape(H1,4,1), reshape(H2,4,1)];
%                     b = reshape(eye(2), 4, 1);
%                 end
%                 if (length(Lt) == 3)
%                     H1 = Lt(1,:)' * Lt(1,:);
%                     H2 = Lt(2,:)' * Lt(2,:);
%                     H3 = Lt(3,:)' * Lt(3,:);
%                     A = [reshape(H1,9,1), reshape(H2,9,1), reshape(H3,9,1)];
%                     b = reshape(eye(3), 9, 1);
%                 end
%                 AMatrix( (((ii-1)*n_elements+1):ii*n_elements), :) = A;
%                 bVector( (((ii-1)*n_elements+1):ii*n_elements) ) = b;
%             end
%             Lambda = (AMatrix' * AMatrix) \ (AMatrix' * bVector);
%             if (nnz(Lambda < 0))
%                 fprintf(2, "Error@ There exist negative components in Lambda\n");
%                 Lambda = abs(Lambda);
%             end
%             this.mSqrtLambda = sqrt(Lambda);
%             if (sum(this.EuclideanTransformationSign < 0) > sum(this.EuclideanTransformationSign > 0))
%                 this.mShape(1,:) = - this.mShape(1,:); % flip the sign of the first row.
%                 this.setWeight(); % weights are updated with new shape
%                 this.EuclideanTransformationSign = - this.EuclideanTransformationSign;
%             end
%         end
        
        
        
        % min  \sum_t || R Pt + t - SqrtLambda * Qt ||_F^2
        function [norm_update, cost] = PoseExtractionByNonlinearLeastSquaresOneIteration (this)
            if this.dim == 2
                rot_dim = 1;
                n_elements = 5 * this.numPointClouds + 4;
            end
            if this.dim == 3
                rot_dim = 3;
                n_elements = 9 * (3*this.numPointClouds + 1);
            end
            this.TmpA_ijv.i_vec = zeros(1, n_elements);
            this.TmpA_ijv.j_vec = zeros(1, n_elements);
            this.TmpA_ijv.v_vec = zeros(1, n_elements);
            vecrhs = zeros(rot_dim*this.numPointClouds+this.dim, 1);
            BB = zeros(this.dim, this.dim);
            bb = zeros(this.dim, 1);
            buffer_cur_pos = 1;
            for ii = 1 : this.numPointClouds
                AA = zeros(rot_dim, rot_dim);
                AB = zeros(rot_dim, this.dim);
                RHS = zeros(rot_dim, 1);
                Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));
                Qt = this.mShape(:, logical(this.PointClouds(ii).Vis));
                %Qt = full(this.transformPoints(Pt, ii));
                %Qt = full(this.mWeights{ii}' * this.FeatureClouds(ii).FeatureMatr(:, logical(this.PointClouds(ii).Vis)));
                meanQt = mean(Qt, 2); meanPt = mean(Pt, 2);
                R = this.mPoses{ii}(:, 1:this.dim);
                for kk = 1 : size(Pt, 2)
                    zcPtk = Pt(:, kk) - meanPt;
                    zcQtk = Qt(:, kk) - meanQt;
                    vrhsk = R * zcPtk - diag(this.mSqrtLambda) * zcQtk;
                    if this.dim == 3
                        Ak = SO3.Hat(R * zcPtk);
                    end
                    if this.dim == 2
                        Ak = [0, 1; -1, 0] * R * zcPtk;
                    end
                    Bk = diag(zcQtk);
                    AA = AA + Ak' * Ak;
                    AB = AB + Ak' * Bk;
                    BB = BB + Bk' * Bk;
                    RHS = RHS + Ak' * vrhsk;
                    bb = bb + Bk' * vrhsk;
                end
                vecrhs(((ii-1)*rot_dim+1) : ii*rot_dim) = RHS;
                buffer_cur_pos = this.putMatrixBlockToStorage(AA, (ii-1)*rot_dim+1, (ii-1)*rot_dim+1, buffer_cur_pos);
                buffer_cur_pos = this.putMatrixBlockToStorage(AB, (ii-1)*rot_dim+1, this.numPointClouds*rot_dim+1, buffer_cur_pos);
                buffer_cur_pos = this.putMatrixBlockToStorage(AB', this.numPointClouds*rot_dim+1, (ii-1)*rot_dim+1, buffer_cur_pos);
            end
            vecrhs(this.numPointClouds*rot_dim+1 : end) = bb;
            buffer_cur_pos = this.putMatrixBlockToStorage(BB, this.numPointClouds*rot_dim+1, this.numPointClouds*rot_dim+1, buffer_cur_pos);     
            update = sparse(this.TmpA_ijv.i_vec, this.TmpA_ijv.j_vec, this.TmpA_ijv.v_vec) \ vecrhs;
            % update rotation, and SqrtLambda
            for ii = 1 : this.numPointClouds
                w = update(((ii-1)*rot_dim+1) : ii*rot_dim);
                if this.dim == 3
                    R = SO3.Exp(w) * this.mPoses{ii}(:, 1:this.dim);
                end
                if this.dim == 2
                    R = [cos(w), -sin(w); sin(w), cos(w)] * this.mPoses{ii}(:, 1:this.dim);
                end
                this.mPoses{ii}(:, 1:this.dim) = R;
            end
            this.mSqrtLambda = this.mSqrtLambda + update(this.numPointClouds*rot_dim+1 : end);
            % detection reflection if sign of this.mSqrtLambda is negative
            reflection_sign = sign(this.mSqrtLambda);
            if sum(reflection_sign < 0) > 0
                this.mSqrtLambda = reflection_sign .* this.mSqrtLambda;
                this.mShape = diag(reflection_sign) * this.mShape;
                this.setWeight(); % update weight parameters
            end
            norm_update = norm(update, 'fro');
            % compute the new cost, and update translation
            cost = 0;
            num_vis = 0;
            for ii = 1 : this.numPointClouds
                num_vis = num_vis + nnz(this.PointClouds(ii).Vis);
                Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));
                Qt = this.mShape(:, logical(this.PointClouds(ii).Vis));
                %Qt = full(this.transformPoints(Pt, ii));
                %Qt = full(this.mWeights{ii}' * this.FeatureClouds(ii).FeatureMatr(:, logical(this.PointClouds(ii).Vis)));
                meanQt = mean(Qt, 2); meanPt = mean(Pt, 2);
                R = this.mPoses{ii}(:, 1:this.dim);
                t =  - (R * meanPt - diag(this.mSqrtLambda) * meanQt);
                this.mPoses{ii}(:, this.dim+1) = t;  % update translation
                vrhs = R * (Pt - meanPt) - diag(this.mSqrtLambda) * (Qt - meanQt);
                normtmp = norm(vrhs, 'fro');
                cost = cost + normtmp * normtmp;
            end
            cost = sqrt(cost/num_vis);
        end
        
        
        
        function [buffer_cur_pos] = putMatrixBlockToStorage(this, A, row, col, buffer_cur_pos)
            
            [m, n] = size(A);
            
            buffer_end = m * n + buffer_cur_pos - 1;
            
            this.TmpA_ijv.i_vec(buffer_cur_pos : buffer_end) = kron (ones(1,n), row:row+m-1);
            
            this.TmpA_ijv.j_vec(buffer_cur_pos : buffer_end) = kron (col:col+n-1, ones(1,m));
            
            this.TmpA_ijv.v_vec(buffer_cur_pos : buffer_end) = full (reshape (A, 1, m * n));
            
            buffer_cur_pos = buffer_end + 1;
            
        end
        
        
    end
    



    methods(Access = private, Static = true)

        function error_statistics = PointCloudsVariationError (PointCloudsCell, VisVecCell)

            num_pointClouds = length(PointCloudsCell);

            array_num_points = zeros(1, num_pointClouds);
            for ii = 1 : num_pointClouds
                array_num_points(ii) = length(PointCloudsCell{ii});
            end
            if ( min(array_num_points) ~= max(array_num_points) )
                fprintf(2, 'Different number of points in each point-clouds \n');
            else
                num_points = array_num_points(1);
            end

            dim = size(PointCloudsCell{1}, 1);
            ref_point_cloud = zeros(dim, num_points);
            sum_vis = zeros(1, num_points);

            for ii = 1 : num_pointClouds
                ref_point_cloud = ref_point_cloud + PointCloudsCell{ii} .* logical(VisVecCell{ii});
                sum_vis = sum_vis + logical(VisVecCell{ii});
            end

            for jj = 1 : length(sum_vis)
                if (sum_vis(jj) > 0)
                    ref_point_cloud(:, jj) = ref_point_cloud(:, jj) ./ sum_vis(jj);
                end
            end

            error_statistics = 0;
            for ii = 1 : num_pointClouds
                Points = PointCloudsCell{ii};
                vis = logical(VisVecCell{ii});
                diff = norm((Points - ref_point_cloud) .* vis, 'fro');
                error_statistics = error_statistics + diff * diff;
            end
            error_statistics = sqrt( error_statistics/sum(sum_vis) );

        end


    end
    
    
end

