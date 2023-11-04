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


classdef KernelGPA < handle
    
    
    % input properties
    properties (Access = public)
        
        dim = 3
        
        PointClouds = struct ( 'Data', {}, 'Vis', {}, 'TransPrior', {})
        
        numPointClouds = 0
        
        numPoints = 0
        
        centerDataAhead = false
        
        flagJointlyOptimizeLamdaRigidTransformation = true;
        
        smoothParam = 0.05 % param for regularizations of warps
        
        kernelScaleParam = 0;

        quantilePercentParam = 0;

        verbose = true;
        
    end
    
    
    % output properties
    properties (Access = public)
        
        mShape = []
        
        mWeights = {}
        
        mPoses = {}
        
        mRsd = []
        
        mTime = []
        
        mSqrtLambda = []
        
        mReferenceCovariancePrior = []
        
        mAffineTransform = {}

        mKernelParams = [];

        mRegularizationStrength = []; % smoothParam * numPoints
        
    end
    
    
    % internal private properties
    properties (Access = private)
        
        H = {}  %  (I-P) * K * inv(K*(I-P)*K + mu*K)
        
        Q = [];
        
        maxNum = 1000000  % the celing number, should never be reached.
        
        EuclideanTransformationSign = []
        
    end
    
    
    methods(Access = public, Static = true)
        
        hFig = VisualizeCameraAndPointclouds3D (fid, PoseCell, Points3D)
       
        [axfig, error_statistics] = VisualizePointCloudsVariation (htfig, PointCloudsCell)
        
    end
    
    
    
    methods (Access = public)
        
        % class constructor
        function this = KernelGPA ( modelFlag )
            if (nargin > 0)
            end
        end
        
        
        
        function addDataByArray(this, DataCell, VisVecCell)
            for ii = 1 : size(DataCell,2)
                this.addPointCloud (DataCell{ii}, VisVecCell{ii});
            end
        end
        
        
        
        function run (this, flagJointlyOptimizeLamdaRigidTransformation)
            if exist('flagJointlyOptimizeLamdaRigidTransformation', 'var') && ~isempty(flagJointlyOptimizeLamdaRigidTransformation)
                this.flagJointlyOptimizeLamdaRigidTransformation = flagJointlyOptimizeLamdaRigidTransformation;
            end            
            this.RunFullProcedure();            
        end
        
        
        function nonrigid_transform_cost_ref = rsdError (this)
            % The reference space cost.
            nonrigid_transform_cost_ref = 0;
            num_vis = 0;
            for ii = 1 : this.numPointClouds
                num_vis = num_vis + nnz(this.PointClouds(ii).Vis);
                pValArray = this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis));
                pValArray = this.transformPoints(pValArray, ii) - this.mShape(:, logical(this.PointClouds(ii).Vis));
                normtmp = norm(pValArray, 'fro');
                nonrigid_transform_cost_ref = nonrigid_transform_cost_ref + normtmp * normtmp;
            end
            nonrigid_transform_cost_ref = sqrt(nonrigid_transform_cost_ref/num_vis);
            this.mRsd = nonrigid_transform_cost_ref;
        end
        
        
        function pImageArray = transformPoints (this, pValArray, cloudIndex )
            pImageArray = pValArray + this.PointClouds(cloudIndex).TransPrior;
            qK = this.CalculateKernelVector(cloudIndex, pImageArray);
            vis = logical(this.PointClouds(cloudIndex).Vis);
            A_t = this.mAffineTransform{cloudIndex};
            A = A_t(1:this.dim, 1: this.dim);
            t = A_t(1:this.dim, this.dim+1);
            pImageArray = (A * pImageArray + t) + this.mShape(:, vis) * this.H{cloudIndex} * qK;
        end


        function pImageArray = affineTransformPoints (this, pValArray, cloudIndex )
            pImageArray = pValArray + this.PointClouds(cloudIndex).TransPrior;
            A_t = this.mAffineTransform{cloudIndex};
            A = A_t(1:this.dim, 1: this.dim);
            t = A_t(1:this.dim, this.dim+1);
            pImageArray = (A * pImageArray + t);
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


        function K = getKernelMatrix(this, cloudIndex, testkernelScaleParam)
            Pts = this.PointClouds(cloudIndex).Data(:, logical(this.PointClouds(cloudIndex).Vis));
            % downsample point-cloud to 5%
            pc = pointCloud(Pts');
            pc = pcdownsample(pc, 'random', 0.05);
            Pts = pc.Location';
            % compute pairwise distances
            K = zeros(size(Pts, 2));
            for kk = 1 : size(Pts, 2)
                K(:, kk) = sum((Pts - Pts(:, kk)).*(Pts - Pts(:, kk)),  1)'; % sum row-wise
            end
            K = sqrt(K);
            mask = triu(true(size(K)),1); % elements above the diagonal
            sigma = median(K(mask)) * testkernelScaleParam;

            % construct K for verification
            K = zeros(size(Pts, 2));
            for kk = 1 : size(Pts, 2)
                K(:, kk) = this.GaussianKernel(Pts, Pts(:, kk), sigma);
            end

        end



        function arrayTrainingStatistics = OptimizeHyperParams(this, smoothParams,  kernelScaleParams, Nfold)
            if ~exist('smoothParams', 'var') || (exist('smoothParams', 'var') && isempty(smoothParams))
                smoothParams = [0.001, 0.01, 0.1, 1];
            end
            if ~exist('kernelScaleParams', 'var') || (exist('kernelScaleParams', 'var') && isempty(kernelScaleParams))
                kernelScaleParams = [0.1, 0.25, 0.5,  0.75,  1.0,  1.25,  1.5, 1.75, 2];
            end
            if ~exist('Nfold', 'var') || (exist('Nfold', 'var') && isempty(Nfold))
                Nfold = 5;
            end
            arrayTrainingStatistics.SmoothParam = zeros(length(smoothParams),  length(kernelScaleParams));
            arrayTrainingStatistics.KernelScaleParam = zeros(length(smoothParams),  length(kernelScaleParams));
            arrayTrainingStatistics.CVE = zeros(length(smoothParams),  length(kernelScaleParams));
            arrayTrainingStatistics.RMSE_ref = zeros(length(smoothParams),  length(kernelScaleParams));
            for ii = 1 : length(smoothParams)
                for jj = 1 : length(kernelScaleParams)
                    this.smoothParam = smoothParams(ii);
                    this.kernelScaleParam = kernelScaleParams(jj);
                    this.RunFullProcedure();
                    RMSE_ref = this.rsdError();
                    cve = this.run_N_Fold_CrossValidation (Nfold);
                    arrayTrainingStatistics.SmoothParam(ii, jj) = smoothParams(ii);
                    arrayTrainingStatistics.KernelScaleParam(ii, jj) = kernelScaleParams(jj);
                    arrayTrainingStatistics.RMSE_ref(ii, jj) = RMSE_ref;
                    arrayTrainingStatistics.CVE(ii, jj) = cve;                    
                end
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
        end        

    end
    
    
    
    
    
    
    % the overall procedure to solve one problem instance
    methods (Access = private)
        


        function setTuningParams (this)
             medianPairwiseDistances = zeros(1, this.numPointClouds);
             meanPairwiseDistances = zeros(1, this.numPointClouds);
             quantileDistances = zeros(1, this.numPointClouds);
             % set kernel parameters
             for ii = 1 : this.numPointClouds
                 Pts = this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis));
                 % downsample point-cloud to 5%
                 if(0)
                     pc = pointCloud(Pts');
                     pc = pcdownsample(pc, 'random', 0.05);
                     Pts = pc.Location';
                 end
                 % compute pairwise distances
                 K = zeros(size(Pts, 2));
                 for kk = 1 : size(Pts, 2)
                     K(:, kk) = sum((Pts - Pts(:, kk)).*(Pts - Pts(:, kk)),  1)'; % sum row-wise
                 end
                 K = sqrt(K);
                 mask = triu(true(size(K)),1); % elements above the diagonal
                 medianPairwiseDistances(ii) = median(K(mask));
                 meanPairwiseDistances(ii) = mean(K(mask));
                 quantileDistances(ii) = quantile(K(mask), this.quantilePercentParam);
             end
             %median_mean_distances = [medianPairwiseDistances; meanPairwiseDistances]
             for ii = 1 : this.numPointClouds
                 if this.kernelScaleParam > 0
                    this.mKernelParams(ii) = meanPairwiseDistances(ii) * this.kernelScaleParam;
                 else
                    this.mKernelParams(ii) = quantileDistances(ii);
                 end
                 this.mRegularizationStrength(ii) = this.smoothParam; % * nnz(this.PointClouds(ii).Vis);
             end
        end



        function RunFullProcedure (this)
             tic;

             % reset data structures
                this.mShape = [];
                this.mWeights = {};
                this.mPoses = {};
                this.mRsd = [];
                this.mTime = [];
                this.mSqrtLambda = [];
                this.mReferenceCovariancePrior = [];
                this.mAffineTransform = {};
                this.H = {};  %  (I-P) * K * inv(K*(I-P)*K + mu*K)
                this.Q = [];
                this.TmpA_ijv = struct('i_vec', [], 'j_vec', [], 'v_vec', []);

             if (length(this.mKernelParams)~=this.numPointClouds) || (length(this.mRegularizationStrength)~=this.numPointClouds)
                 this.setTuningParams();
             end

                       
             if (~this.flagJointlyOptimizeLamdaRigidTransformation)
                this.calReferenceCovariancePrior('ascend'); % ascending order
             else
                this.mReferenceCovariancePrior = ones(this.dim, 1);
                transPrior = zeros(this.dim, 1);
                for ii = 1 : this.numPointClouds
                    this.PointClouds(ii).TransPrior = transPrior * logical(this.centerDataAhead);
                end
            end
            
            % construct Q to factorize
            this.Q = zeros(this.numPoints, this.numPoints);
            for ii = 1 : this.numPointClouds
                Kt = this.CaluclateKernelMatrix(ii);
                proj_I_minus_Pt = this.CalculateOrthogonalProjectionMatrix (ii);
                this.H{ii} = proj_I_minus_Pt * inv(Kt * proj_I_minus_Pt + this.mRegularizationStrength(ii) * eye(length(Kt)));
                % construct Qt, for each point-cloud
                nsize = size(Kt, 2);
                Qt = - (this.H{ii} * Kt - eye(nsize)) * proj_I_minus_Pt;
                Qt = (Qt + Qt')/2;
                % construct:  Q = sum_t Tau * Qt * Tau'.   With Tau being  truncated visiblity matrix
                % Tau = [... e_k ...].   If k-th point is visible in this point-cloud
                vis = logical(this.PointClouds(ii).Vis);
                this.Q(vis, vis) = this.Q(vis, vis) + Qt;

                if(0)
                % BEGIN TEST EQUALITY
                oldQt = (this.H{ii} * Kt - eye(nsize)) * proj_I_minus_Pt;
                oldQt = oldQt * oldQt' + this.mRegularizationStrength(ii) * this.H{ii} * Kt * this.H{ii}';                
                % END TEST EQUALITY
                %fprintf(2, "%d Negative components in: St \n", nnz(eig((St+St')/2) < 0));                
                %fprintf(2, "%d Negative components in: I - Kt_invSt_Kt \n", nnz(eig((Yt+Yt')/2) < 0));
                %fprintf(2, "%d Negative components in: I - (I-P)Kt_invSt_Kt(I-P) \n", nnz(eig((YYt+YYt')/2) < 0));
                %fprintf(2, "%d Negative components in: Qt = (I-P) - (I-P)Kt_invSt_Kt(I-P) \n", nnz(eig((Qt+Qt')/2) < 0));
                norm0 = norm((this.H{ii}-this.H{ii}'), 'fro');
                norm1 = norm((this.mRegularizationStrength(ii)*this.H{ii}-Qt), 'fro');
                norm2 = norm((this.mRegularizationStrength(ii)*this.H{ii}-oldQt), 'fro');
                norm3 = norm((Qt-oldQt), 'fro');
                if norm0 > 1e-12 || norm1 > 1e-12 || norm2 > 1e-12 || norm3 > 1e-12
                fprintf(2, 'norm(H - H^t) = %.25f \n', norm0);
                fprintf(2, 'norm(mu*H - Qt) = %.25f \n', norm1);
                fprintf(2, 'norm(mu*H - oldQt) = %.25f \n', norm2);
                fprintf(2, 'norm(Qt - oldQt) = %.25f \n', norm3);
                end
%                eig(Qt)
                end

            end
            if norm(this.Q * ones(length(this.Q), 1)) > 1e-6
                fprintf(2, 'Fatal Error. Q*1=0 is violated ||Q*1|| = %.15f \n',   norm(this.Q * ones(length(this.Q), 1)) );
                %return;
            end
            
            % Q * 1 = 0.
            this.Q = this.Q + (this.numPointClouds); % soft regularization, to avoid eigenvector 1.
            [EV, eigval, ~] = this.byEigen (this.Q, this.dim, 'smallest'); % descending order
            this.mShape =  EV';

            if ~this.flagJointlyOptimizeLamdaRigidTransformation
                this.eliminateReflection();
            end

            % Initialize the affine part of the transformation model
            for ii = 1 : this.numPointClouds
                Kt = this.CaluclateKernelMatrix(ii);
                nsize = length(Kt);
                Pts = [this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis));  ones(1, nsize)];
                this.mAffineTransform{ii} = - this.mShape(:, logical(this.PointClouds(ii).Vis)) * (this.H{ii} * Kt - eye(nsize)) * Pts' * pinv(Pts * Pts');
            end
            
            
            if this.flagJointlyOptimizeLamdaRigidTransformation

%                Lambda = this.LambdaFromRigidAffine();
%                fprintf(1, "sqrtLambda = ( %f ", this.mSqrtLambda); fprintf(1, ")\n");

                this.InitializeLambda();
                lambda_rigid_transform_cost = this.InitializePoses();
                if (this.verbose)
                fprintf(1, 'iter[0]   cost = %f  sqrtLambda = (', lambda_rigid_transform_cost);  fprintf(1, ' %f ', this.mSqrtLambda); fprintf(1, ')\n');
                end
                for ii = 1 : 100
                    [norm_update, lambda_rigid_transform_cost] = this.PoseExtractionByNonlinearLeastSquaresOneIteration();
                    if (this.verbose)
                    fprintf(1, 'iter[%d]  cost = %f  ||update|| = %f  sqrtLambda = (', ii, lambda_rigid_transform_cost, norm_update);  fprintf(1, ' %f ', this.mSqrtLambda); fprintf(1, ')\n');
                    end
                    if norm_update < 1e-7
                        break;
                    end
                end
                this.mReferenceCovariancePrior = this.mSqrtLambda .* this.mSqrtLambda;
                % scale the reference shape and affine transformation with the estimated SqrtLambda
                this.mShape = diag(this.mSqrtLambda) * this.mShape;
                for ii = 1 : this.numPointClouds
                    this.mAffineTransform{ii} = diag(this.mSqrtLambda) * this.mAffineTransform{ii};
                end
            else
                this.mSqrtLambda = sqrt(this.mReferenceCovariancePrior);
                % scale the reference shape and affine transformation with the estimated SqrtLambda
                this.mShape = diag(this.mSqrtLambda) * this.mShape;
                for ii = 1 : this.numPointClouds
                    this.mAffineTransform{ii} = diag(this.mSqrtLambda) * this.mAffineTransform{ii};
                end
                this.ConstructOptimalRigidTransform();
            end
            
            %this.RoundPoseFromAffine(); % Give Lambda, solve min || A - R ||_F^2
            %this.ConstructOptimalRigidTransform(); % Give Lambda, solve || R Pt + t - y_t(P_t) ||
            %this.DecomposePoseAndDeformation();

            this.mTime(1) = toc;
            
        end


        function [transform_points_cell, transform_vis_cell, refShape] = runTrainingTestData(this, trainingPointsFlag, testingPointsFlag)
            AllPointClouds = this.PointClouds;
            AllNumPoints = this.numPoints;
            this.PointClouds = struct ( 'Data', {}, 'Vis', {}, 'TransPrior', {});
            this.numPointClouds = 0;
            this.numPoints = 0;
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
        
        
    end
    
    
    
    
    

    



    methods (Access = private)


        function DecomposePoseAndDeformation (this)

            % set kernel parameters
            Pts = this.mShape;
            % downsample point-cloud to 5%
            pc = pointCloud(Pts');
            pc = pcdownsample(pc, 'random', 0.05);
            Pts = pc.Location';
            % compute pairwise distances
            K = zeros(size(Pts, 2));
            for kk = 1 : size(Pts, 2)
                K(:, kk) = sum((Pts - Pts(:, kk)).*(Pts - Pts(:, kk)),  1)'; % sum row-wise
            end
            K = sqrt(K);
            mask = triu(true(size(K)),1); % elements above the diagonal

            msigma = median(K(mask)) * this.kernelScaleParam;
            mu = 2;


            % construct the kernel matrix
            mK = zeros(size(this.mShape, 2));            
            for kk = 1 : size(this.mShape, 2)
                mK(:, kk) = this.GaussianKernel(this.mShape, this.mShape(:, kk), msigma);
            end

            for ii = 1 : this.numPointClouds

                Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));

                Qt = this.mShape(:, logical(this.PointClouds(ii).Vis)); 

                vis = logical(this.PointClouds(ii).Vis);

                nsize = nnz (vis);

                % min || R * Pt + t * 1' - Qt - Wt * mK ||_F^2  ---> (R, t)

                invmK = inv( mK(vis, vis) + mu * eye(nsize) );
                Nt = eye(nsize) - invmK;

                % min || (R * Pt + t * 1' - Qt) * Nt ||_F^2  ---> (R, t)

                Nt1 = Nt' * ones(nsize, 1);
                ProjNt1 = eye(nsize) - (Nt1 * Nt1.')./(Nt1.' * Nt1);



                % min || (R * Pt * Nt - Qt * Nt) * ProjNt1 ||_F^2
                PtBar = Pt * Nt * ProjNt1;
                QtBar = Qt * Nt * ProjNt1;

                % min ||R * PtBar - QtBar ||_F^2
                R = (PtBar - mean(PtBar, 2)) * (QtBar - mean(QtBar, 2))';
                [U, ~, V] = svd(R); det_VU = det(V*U'); % this has the same sign as det(R). If det(R) > 0, then det(V'U) > 0
                diagS = eye(length(R));  diagS(end, end) = det_VU;
                % the optimal rotation
                R = V * diagS * U';
                % the optimal translation
                t = - (R * Pt * Nt - Qt * Nt) * (Nt1 ./ (Nt1.' * Nt1));
                
                % the optimal pose
                this.mPoses{ii} = [R, t];

                % the deformation field: Qt + Wt * mK
                Wt = ( (R * Pt + t) - Qt ) * invmK;

            end

        end



    end








    properties (Access = private)
        
        TmpA_ijv = struct('i_vec', [], 'j_vec', [], 'v_vec', []);
        
    end
    


    % methods used to compute Lambda and poses
    methods (Access = private)
        



        function Lambda = LambdaFromRigidAffine (this)
            this.EuclideanTransformationSign = zeros(1, this.numPointClouds);
            RgInvLambdaRg = zeros(this.dim ,this.dim);
            for ii = 1 : this.numPointClouds
                A_t = this.mAffineTransform{ii};
                Lt = A_t(1:this.dim, 1: this.dim);
                this.EuclideanTransformationSign(ii) = sign(det(Lt));  % the sign of affine transformation
                %  || L * L' - Inv(Lambda) ||_F = || diag( L * L' - Inv(Lambda) ) ||_2
                RgInvLambdaRg = RgInvLambdaRg + Lt * Lt';       
            end
            
            [Rg, diagInvLambda] = eig(RgInvLambdaRg ./ this.numPointClouds);          
            inv_lambda = diag(diagInvLambda);
            Lambda = 1 ./ inv_lambda;
            % arrange elements of Lambda in the ascending order
            [Lambda, index] = sort(Lambda, 'ascend');
            Rg = Rg(:, index);

            %testouput_lambda = Lambda
            if (nnz(Lambda < 0))
                fprintf(2, 'Error@ There exist negative components in Lambda\n');
                Lambda = abs(Lambda);
            end
            this.mSqrtLambda = sqrt(Lambda);
            if (sum(this.EuclideanTransformationSign < 0) > sum(this.EuclideanTransformationSign > 0))
                this.mShape(1,:) = - this.mShape(1,:); % flip the sign of the first row.
                this.mPoses{ii}(1,:) = - this.mPoses{ii}(1,:);
                this.EuclideanTransformationSign = - this.EuclideanTransformationSign;
            end
        end



        % rounding affine A to the nearest rotation 
        % min || A - R ||_F^2
        function RoundPoseFromAffine (this)
            for ii = 1 : this.numPointClouds
                A_t = this.mAffineTransform{ii};
                A = A_t(1:this.dim, 1: this.dim);
                t = A_t(1:this.dim, this.dim+1);
                [U, S, V] = svd(A);
                if this.dim == 2
                   R = U * diag([1, det(U*V')]) * V';
                end
                if this.dim == 3
                   R = U * diag([1, 1, det(U*V')]) * V';
                end
                this.mPoses{ii} = [R, t];
            end
        end



        % || R Pt + t - y_t(P_t) ||
        function ConstructOptimalRigidTransform (this)
            for ii = 1 : this.numPointClouds
                Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));
                TransPt = full(this.transformPoints(Pt, ii));
                %TransPt = full(this.affineTransformPoints(Pt, ii));
                meanTransPt = mean(TransPt, 2);  meanPt = mean(Pt, 2);
                R = (Pt - meanPt) * (TransPt - meanTransPt)';
                [U, ~, V] = svd(R); det_VU = det(V*U'); % this has the same sign as det(R). If det(R) > 0, then det(V'U) > 0
                diagS = eye(this.dim);  diagS(end, end) = det_VU;
                R = V * diagS * U';
                t = - (R * meanPt - meanTransPt);
                this.mPoses{ii} = [R, t];
            end
        end





        function cost = InitializePoses (this)
            if ((abs(sum(this.EuclideanTransformationSign)) ~= this.numPointClouds)  && this.verbose)
                fprintf(2, 'Consistent optimal rigid transformation to the reference point-cloud DO NOT exist. \n');
                sign_of_euclidean_transformations = this.EuclideanTransformationSign
            end
            for ii = 1 : this.numPointClouds
                Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));
                % Qt = this.mShape(:, logical(this.PointClouds(ii).Vis));
                Qt = full(this.transformPoints(Pt, ii));
                %Qt = full(this.affineTransformPoints(Pt, ii));
                meanQt = mean(Qt, 2);
                meanPt = mean(Pt, 2);
                R = (Pt - meanPt) * (Qt - meanQt)' * diag(this.mSqrtLambda);
                [U, ~, V] = svd(R); det_VU = det(V*U'); % this has the same sign as det(R). If det(R) > 0, then det(V'U) > 0
                diagS = eye(this.dim);  diagS(end, end) = det_VU;
                R = V * diagS * U';
%                 if (det_VU < 0 && this.verbose)
%                     fprintf(2, 'Warning @Det(VU) < 0 in orthogonal Procrustes \n');
%                 end
                t = - (R * meanPt - diag(this.mSqrtLambda) * meanQt);
                this.mPoses{ii} = [R, t];
            end
            % compute the initial cost
            cost = 0;
            num_vis = 0;
            for ii = 1 : this.numPointClouds
                num_vis = num_vis + nnz(this.PointClouds(ii).Vis);
                Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));
                Qt = full(this.transformPoints(Pt, ii));
                R = this.mPoses{ii}(:, 1:this.dim);
                t = this.mPoses{ii}(:, this.dim+1);
                vrhs = (R * Pt + t) - diag(this.mSqrtLambda) * Qt;
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
%                 % Qt = this.mShape(:, logical(this.PointClouds(ii).Vis));
%                 Qt = full(this.transformPoints(Pt, ii));
%                 zcQt = Qt - mean(Qt, 2);
%                 zcPt = Pt - mean(Pt, 2);
%                 Lt = zcQt * zcPt' * pinv(zcPt * zcPt');
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
                %Qt = full(this.affineTransformPoints(Pt, ii));
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
                for ikk = 1 : this.numPointClouds
                    this.mAffineTransform{ikk} = diag(reflection_sign) * this.mAffineTransform{ikk};
                end
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
    
    
    
    
    
    % kernel functions
    methods (Access = private, Static = true)
      
        function val = GaussianKernel (Points, qpt, param)
            cc = -2 * param * param;
            val = sum((Points - qpt).*(Points - qpt),  1)' ./ cc;
            val = exp(val);
        end
        
        
        function val = PolynomialKernel (pt1, pt2, param)
            val = pt1' * pt2;
            val = val^param;
        end
        
    end
    

    % methods used to constructed internal matrices
    methods (Access = private)
        
        
        function qK = CalculateKernelVector (this, cloudIndex, query_pt)
            Pts = this.PointClouds(cloudIndex).Data(:, logical(this.PointClouds(cloudIndex).Vis));
            qK = zeros(size(Pts, 2), size(query_pt, 2));
            sigma = this.mKernelParams(cloudIndex);
            for kk = 1 : size(query_pt, 2)
                qK(:, kk) = this.GaussianKernel(Pts, query_pt(:, kk), sigma);
            end
        end
        
        
        
        function K = CaluclateKernelMatrix (this, cloudIndex)
            Pts = this.PointClouds(cloudIndex).Data(:, logical(this.PointClouds(cloudIndex).Vis));
            K = zeros(size(Pts, 2));
            sigma = this.mKernelParams(cloudIndex);
            for kk = 1 : size(Pts, 2)
                K(:, kk) = this.GaussianKernel(Pts, Pts(:, kk), sigma);
            end
        end
        
        
        
        function proj_I_minus_Pt = CalculateOrthogonalProjectionMatrix (this, cloudIndex)
            Pts = this.PointClouds(cloudIndex).Data(:, logical(this.PointClouds(cloudIndex).Vis));
            nsize = size(Pts, 2);
            Pts = [Pts; ones(1, nsize)];
            proj_I_minus_Pt = eye(nsize) - Pts' * pinv(Pts * Pts') * Pts;
        end
        
    end
    
    
    

    methods (Access = private, Static = true)
        
        function [error_statistics,  normalizer_const] = PointCloudsVariationError (PointCloudsCell,  VisVecCell)
            error_statistics = 0;
            num_pointClouds = length(PointCloudsCell);
            num_points = size(PointCloudsCell{1}, 2);
            % compute mean
            ref_point_cloud = 0 * PointCloudsCell{1};
            ref_num_vis = 0 * VisVecCell{1};
            for ii = 1 : num_pointClouds
                ref_point_cloud = ref_point_cloud + PointCloudsCell{ii} .* VisVecCell{ii};
                ref_num_vis = ref_num_vis + VisVecCell{ii};
            end
            for jj = 1 : num_points
                if (ref_num_vis(jj) > 0)
                    ref_point_cloud(:, jj) = ref_point_cloud(:, jj) ./ ref_num_vis(jj);
                end
            end
            % compute error
            for ii = 1 : num_pointClouds
                Diff = (PointCloudsCell{ii} - ref_point_cloud) .* VisVecCell{ii};
                diff = norm(Diff, 'fro');
                error_statistics = error_statistics + diff * diff;
            end
            normalizer_const = sum(ref_num_vis);
            error_statistics = sqrt(error_statistics/normalizer_const);
        end
        
    end




    % methods used to calcualte hyper parameters
    methods (Access = private)

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
                Pt = this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)); normtmp = norm((R * Pt + t) - this.transformPoints(Pt, ii), 'fro');
                %Diff = (R * this.PointClouds(ii).Data + t) - this.mShape; normtmp = norm( Diff .* this.PointClouds(ii).Vis, 'fro');
                cost = cost + normtmp * normtmp;
            end
            cost = sqrt(cost/num_vis);
        end





    end
    
    
    
    
    % ------------------------- Methods from IJCV 2022 paper ------------------
    % Bai, F., Bartoli, A. Procrustes Analysis with Deformations: A Closed-Form Solution by Eigenvalue Decomposition. Int J Comput Vis 130, 567–593 (2022). 
    % https://doi.org/10.1007/s11263-021-01571-8
    methods (Access = private)
        
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
                        
                        sumD = sumD + (sRij * this.PointClouds(jj).Data.*this.PointClouds(jj).Vis + tij);
                        
                        if (det(sRij) < 0 && this.verbose)
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
            
            this.mReferenceCovariancePrior = v .* v;
            
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
    
    
    
    
    
    
    % static methods of solvers
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


    methods (Access = public)


        function clear (this)

%             this.dim = 3

            this.PointClouds = struct ( 'Data', {}, 'Vis', {}, 'TransPrior', {})
% 
            this.numPointClouds = 0

%             this.numPoints = 0

%             this.centerDataAhead = false

%             this.flagJointlyOptimizeLamdaRigidTransformation = true;

%             this.smoothParam = 0.05 % param for regularizations of warps

%             this.kernelScaleParam = 1.0;
% 
%             this.verbose = true;

            this.mShape = []

            this.mWeights = {}

            this.mPoses = {}

            this.mRsd = []

            this.mTime = []

            this.mSqrtLambda = []

            this.mReferenceCovariancePrior = []

            this.mAffineTransform = {}

            this.mKernelParams = [];

            this.mRegularizationStrength = []; % smoothParam * numPoints

            this.H = {}  %  (I-P) * K * inv(K*(I-P)*K + mu*K)

            this.Q = [];

            this.maxNum = 1000000  % the celing number, should never be reached.

            this.EuclideanTransformationSign = []

            this.TmpA_ijv = struct('i_vec', [], 'j_vec', [], 'v_vec', []);

        end

    end
    
    
    
end

