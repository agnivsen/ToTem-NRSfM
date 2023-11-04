% Copyright 2022 Fang Bai <fang dot bai at yahoo dot com>
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


classdef RigidGPA < handle
    
    
    % input data
    properties (Access = public)
        
        dim = 3
        
        PointClouds = struct ( 'Data', {}, 'Vis', {},  'PtID', {})
        
        numPointClouds = 0
        
        numPoints = 0
        
    end
    
    
    % output properties
    properties (Access = public)
        
        mShape = []
        
        mPoses = {}
        
        mRsd = []
        
        time = []
        
        dimFeatureSpace = -1
        
        verbose = true
        
        flagEliminateReferenceShapeFirst = false;
        
    end
    
    
    
    
    
    
    
    
    methods (Access = public)
        
        % class constructor
        function this = RigidGPA ()
            
            if (nargin > 0)
                
            end
            
        end
        
        
        
        function addDataByArray(this, DataCell, VisVecCell)
            
            for ii = 1 : size(DataCell,2)
                
                this.addPointCloud (DataCell{ii}, VisVecCell{ii});
                
            end
            
        end
        
        
        % interface function to template experiments.
        function run (this, smoothParam, flagJointlyOptimizeLamdaRigidTransformation)
            this.optimize();
        end
        

        % interface function to template experiments.        
        function [cve, predictedPointClouds] = run_N_Fold_CrossValidation (this, N)
            this.optimize();
            gaugeShape = this.mShape;
            predictedPointClouds = cell(1, this.numPointClouds);
            for ii = 1 : this.numPointClouds
                predictedPointClouds{ii} = 0 * this.PointClouds(ii).Data;
            end
            N = min(this.numPoints, N);
            group_size = this.numPoints/N;
            group_start_pos = [round(1 : group_size : this.numPoints),  this.numPoints+1];
            for ii = 1 : N
                lb = group_start_pos(ii);
                ub = group_start_pos(ii+1)-1;
                testingPointsFlag = false(1, this.numPoints);
                testingPointsFlag(lb:ub) = true;
                trainingPointsFlag = ~logical(testingPointsFlag);
                % compute the transformed test points
                [transform_points_cell, transform_vis_cell, refShape] = runTrainingTestData(this, trainingPointsFlag, testingPointsFlag);
                % correct gauge and fill the predicted points
                 [sR, t] = this.SimilarityProcrustes (refShape, gaugeShape(:, trainingPointsFlag));
                for kk = 1 : this.numPointClouds
                    predictedPointClouds{kk}(:, testingPointsFlag) = (sR * transform_points_cell{kk} + t) .* transform_vis_cell{kk};
                end
            end
            % CVE is computed using the discrepancy between the predicted-point-clouds and the global reference shape
            normalizer = 0; cve = 0;
            for ii = 1 : this.numPointClouds
                vis =  this.PointClouds(ii).Vis;
                pointsDeviation = (predictedPointClouds{ii} - gaugeShape) .* vis;
                sqrtError =  norm(pointsDeviation, 'fro');
                cve = cve + sqrtError * sqrtError;
                normalizer = normalizer + nnz(vis);
            end
            cve = sqrt(cve/normalizer);
            % reset the estimate to the status of using all the points.
            this.optimize();
        end


        
        function optimize (this, lambda, maxIters, errorTol, initPosesCell)
            tic;
            
            if ~exist('lambda', 'var') || (exist('lambda', 'var') && isempty(lambda))
                lambda = 0;
            end
            if ~exist('maxIters', 'var') || (exist('maxIters', 'var') && isempty(maxIters))
                maxIters = 50;
            end
            if ~exist('errorTol', 'var') || (exist('errorTol', 'var') && isempty(errorTol))
                errorTol = 1e-5;
            end
            if  exist('initPosesCell', 'var') && (length(initPosesCell) == this.numPointClouds)
                this.mPoses = initPosesCell;
            else
                % initialize the poses by affine GPA
                this.initalizePoseByAffineGPA ();
            end
            if lambda > 0
                fprintf(1, '\nAlgorithm: Levenberg-Marquardt \n');
            else
                fprintf(1, '\nAlgorithm: Gauss-Newton \n');
                lambda = 0;
            end
            fprintf(1, 'Maximum Iterations: %d \n', maxIters);
            fprintf(1, 'Error Tolerence: %f \n', errorTol);
            this.InitializeVariables ();
            this.UpdateReferenceShape();
            rsd = this.rsdError();
            if (this.verbose)
            fprintf(2, '[%d]  fcost  = %f   \n', 0, rsd);
            end
            for ii = 1 : maxIters
                if this.flagEliminateReferenceShapeFirst
                    % this will eliminate reference shape first
                    if (lambda > 0)
                        this.RunLevenbergMarquardtOneIteration(lambda)   % This use LevenbergMarquardt
                    end
                    if (lambda ==0)
                        this.RunGaussNewtonOneIteration();  % This use GaussNewton
                    end
                    this.UpdateReferenceShape();
                    [flcost, nlcost] = this.fLinearCost (lambda);
                else
                    [flcost, nlcost] = this.SparseNonlinearLeastSquaresOneIteration (lambda);
                end
                norm_update = norm(this.UpdateX, 'fro');
                rsd = this.rsdError();                
                if (this.verbose)
                fprintf(1, '[%d]  fcost  = %f   fcost_prev = %f   flcost = %f  ||dx|| = %f   lambda = %f \n', ii, rsd, nlcost, flcost, norm_update, lambda);
                end
                if (norm_update < errorTol)
                    break;
                end
                if (lambda > 0)
                    % for descent, rho > 0;   for ascent, rho < 0
                    rho = (rsd - nlcost) / (flcost - nlcost);  % always flcost < nlcost.  so  flcost - nlcost < 0
                    if rho > 0
                        lambda = lambda/3;
                    else
                        lambda = lambda*2;
                    end
                end
            end
            fprintf('\n');
            
            this.time(1) = toc;
        end
        
        
        
        function pImageArray = transformPoints (this, pValArray, cloudIndex )
            R = this.mPoses{cloudIndex}(:, 1:this.dim);
            t = this.mPoses{cloudIndex}(:, this.dim+1);
            pImageArray = R * pValArray + t;
        end
        
        
        
        function rsd_ref = rsdError(this)
            rsd_ref = 0;
            for ii = 1 : this.numPointClouds
                R = this.mPoses{ii}(:, 1:this.dim);
                t = this.mPoses{ii}(:, this.dim+1);
                Diff = (R * this.PointClouds(ii).Data + t) - this.mShape;
                dr = norm(Diff .* this.PointClouds(ii).Vis, 'fro');
                rsd_ref = rsd_ref + dr * dr;
            end
            v = sum(this.vecSumVis);
            rsd_ref = sqrt(rsd_ref/v);
            this.mRsd = [this.mRsd, rsd_ref];
        end
        
        
        
        function addPointCloud (this, dataMatr, visVec )
            id_vec = 1 : size(dataMatr, 2);
            cnt = this.numPointClouds + 1;
            this.PointClouds(cnt).Data = double(dataMatr);
            this.PointClouds(cnt).Vis = ( sum (abs(dataMatr) < this.maxNum) == size(dataMatr, 1) );
            if (nargin == 3)
                this.PointClouds(cnt).Vis = logical(visVec);
            end
            this.PointClouds(cnt).PtID = id_vec(logical(visVec));
            this.numPointClouds = cnt;
            this.dim = size(dataMatr, 1);
            this.numPoints = size(dataMatr, 2);
        end
        
        
        
        
        
        
        
    end
    
    
    
    % internal private properties
    properties (Access = private)
        
        maxNum = 1000000  % the celing number, should never be reached.
        
    end
    
    
    
    % temporary objects used in reduced least squares procedure
    properties (Access = private)
        
        HessianMatrix = []
        RhdResidual = []
        
        UpdateX = []
        
        vecSumVis = []
        
        vecSumVisInv = []
        
        
        TmpSquareCk = []
        TmpSquaredk = []
        
        TmpPData = []
        TmpqVis = []
        
        TmpAk
        TmpBk
        
        param_dim_rotation = 3
        
        Tmpbb = 0
        
    end
    
    
    % private methods that used to solve the problem
    methods (Access = private)
        
        
        function RunLevenbergMarquardtOneIteration (this, lambda)
            param_offset = this.param_dim_rotation + this.dim;
            this.HessianMatrix = zeros(param_offset*this.numPointClouds);
            this.RhdResidual = zeros (param_offset*this.numPointClouds, 1);
            this.Tmpbb = 0;
            for ii = 1 : this.numPointClouds
                this.ProcessOneReducedSquare(ii);
            end
            for ii = 1 : length(this.HessianMatrix)
                this.HessianMatrix(ii, ii) = this.HessianMatrix(ii, ii) + lambda;
            end
            this.UpdateX = this.HessianMatrix \ this.RhdResidual;
            for ii = 1 : this.numPointClouds
                deltaX = this.UpdateX( (param_offset*(ii-1)+1) : (param_offset*ii) );
                deltaX_rota = deltaX(1:this.param_dim_rotation);
                if this.dim == 2
                    deltaR = [cos(deltaX_rota),  -sin(deltaX_rota);  sin(deltaX_rota),  cos(deltaX_rota)];
                end
                if this.dim == 3
                    deltaR = SO3.Exp(deltaX_rota);
                end
                R = this.mPoses{ii}(:, 1:this.dim);
                t = this.mPoses{ii}(:, this.dim+1);
                this.mPoses{ii} = [deltaR*R,  t+deltaX((this.param_dim_rotation+1) : end)];
            end
        end
        
        
        
        function RunGaussNewtonOneIteration (this)
            param_offset = this.param_dim_rotation + this.dim;
            this.HessianMatrix = zeros(param_offset*this.numPointClouds);
            this.RhdResidual = zeros (param_offset*this.numPointClouds, 1);
            this.Tmpbb = 0;
            for ii = 1 : this.numPointClouds
                this.ProcessOneReducedSquare(ii);
            end
            % this fixes the first pose
            this.UpdateX = this.HessianMatrix((param_offset+1):end, (param_offset+1):end) \ this.RhdResidual((param_offset+1):end);
            % update poses except the first one
            for ii = 2 : this.numPointClouds
                deltaX = this.UpdateX( (param_offset*(ii-2)+1) : (param_offset*(ii-1)) );
                deltaX_rota = deltaX(1:this.param_dim_rotation);
                if this.dim == 2
                    deltaR = [cos(deltaX_rota),  -sin(deltaX_rota);  sin(deltaX_rota),  cos(deltaX_rota)];
                end
                if this.dim == 3
                    deltaR = SO3.Exp(deltaX_rota);
                end
                R = this.mPoses{ii}(:, 1:this.dim);
                t = this.mPoses{ii}(:, this.dim+1);
                this.mPoses{ii} = [deltaR*R,  t+deltaX((this.param_dim_rotation+1) : end)];
            end
        end
        
        
        
        function UpdateReferenceShape (this)
            this.mShape = zeros(this.dim,  this.numPoints);
            for ii = 1 : this.numPointClouds
                R = this.mPoses{ii}(:, 1:this.dim);
                t = this.mPoses{ii}(:, this.dim+1);
                this.mShape = this.mShape + ( (R * this.PointClouds(ii).Data + t) .* this.PointClouds(ii).Vis );
            end
            this.mShape = this.mShape .* this.vecSumVisInv;
        end
        
        
        
        % this is to obtain Ck and dk
        % used as: Ck' * Ck      and      Ck' * dk
        function ProcessOneReducedSquare (this, ii)
            this.TmpSquareCk = zeros(this.dim*this.numPoints, (this.param_dim_rotation+this.dim)*this.numPointClouds);
            this.TmpSquaredk = zeros(this.dim*this.numPoints, 1);
            % vis = inv(Tau_sum) * Tau_i
            vis = this.vecSumVisInv .* this.PointClouds(ii).Vis;
            param_offset = this.param_dim_rotation + this.dim;
            for kk = 1 : this.numPointClouds
                % Tau_k   * inv(Tau_sum) * Tau_i = Tau_k * vis
                this.TmpqVis = - this.PointClouds(kk).Vis .* vis;
                if (kk == ii)
                    this.TmpqVis = this.TmpqVis + this.PointClouds(ii).Vis;
                end
                this.TmpPData = this.PointClouds(kk).Data .* this.TmpqVis;
                % TmpPData and Tmpqvis
                R = this.mPoses{kk}(:, 1:this.dim);
                t = this.mPoses{kk}(:, this.dim+1);
                this.TmpPData = R * this.TmpPData;
                % Rk * Pk + tk * 1^T * Qk
                % Get the linearized coefficients:   Ak,   Bk
                this.GetRotatonTransaltionCoefficient();
                this.TmpSquareCk(:,  ((kk-1)*param_offset+1) : kk*param_offset) = [ this.TmpAk, this.TmpBk ];
                % Get residual added up:       dk = h1+h2+ ... + hn
                this.TmpSquaredk = this.TmpSquaredk + reshape((this.TmpPData + t*this.TmpqVis), this.dim*this.numPoints ,1);
            end
            this.HessianMatrix = this.HessianMatrix + this.TmpSquareCk'*this.TmpSquareCk;
            this.RhdResidual = this.RhdResidual - this.TmpSquareCk'*this.TmpSquaredk;
            this.Tmpbb = this.Tmpbb + this.TmpSquaredk' * this.TmpSquaredk;
        end

                
        
        function GetRotatonTransaltionCoefficient (this)
            % linearize rotation
            if this.dim == 2
                for ii = 1 : length(this.TmpPData)
                    this.TmpAk((ii*this.dim-this.dim+1) : ii*this.dim,  1:this.param_dim_rotation) = [0 , -1; 1, 0] * this.TmpPData(:, ii);
                end
            end
            if this.dim == 3
                for ii = 1 : length(this.TmpPData)
                    this.TmpAk((ii*this.dim-this.dim+1) : ii*this.dim,  1:this.param_dim_rotation) = -SO3.Hat(this.TmpPData(:, ii));
                end
            end
            % linerize translation
            for ii = 1 : length(this.TmpPData)
                this.TmpBk((ii*this.dim-this.dim+1) : ii*this.dim,  1:this.dim) = this.TmpqVis(ii) * eye(this.dim);
            end
        end
        
        
        
        % ||Ax - b|| = x' A'*A x - 2*x'* A'*b + b'*b
        function [flcost, nlcost] = fLinearCost (this, lambda)            
            if (lambda > 0)
                xAAx = this.UpdateX' * (this.HessianMatrix - lambda * speye(length(this.HessianMatrix))) * this.UpdateX;
                xAb = this.UpdateX' * this.RhdResidual;
            end
            if (lambda == 0)
                param_offset = this.param_dim_rotation + this.dim;
                xAAx = this.UpdateX' * this.HessianMatrix((param_offset+1):end, (param_offset+1):end) * this.UpdateX;
                xAb = this.UpdateX' * this.RhdResidual((param_offset+1):end);
            end
            bb = this.Tmpbb;
            flcost = xAAx - 2*xAb + bb;
            v = sum(this.vecSumVis);
            nlcost = sqrt(bb/v);
            flcost = sqrt(flcost/v);
        end
        
        
        
        function InitializeVariables (this)
            if this.dim == 2
                this.param_dim_rotation = 1;
            end
            if this.dim == 3
                this.param_dim_rotation = 3;
            end
            this.vecSumVis = zeros(1, this.numPoints);
            for ii = 1 : this.numPointClouds
                this.vecSumVis = this.vecSumVis + this.PointClouds(ii).Vis;
            end
            this.TmpSquareCk = zeros(this.dim*this.numPoints, (this.param_dim_rotation+this.dim)*this.numPointClouds);
            this.TmpSquaredk = zeros(this.dim*this.numPoints, 1);
            this.TmpPData = zeros(this.dim, this.numPoints);
            this.TmpqVis = zeros(1, this.numPoints);
            this.TmpAk = zeros(this.dim*this.numPoints, this.param_dim_rotation);
            this.TmpBk = zeros(this.dim*this.numPoints, this.dim);
            this.vecSumVisInv = 1 ./ this.vecSumVis;
            this.mRsd = [];
        end
        
        
    end
    
    
    
    % this implements sparse linear squares
    
    
    properties (Access = private)

        SparseHessianMatrix_ijv = struct('i_vec', [], 'j_vec', [], 'v_vec', []);
        
        SparseRhsVector = []
        
        
    end
    
    methods (Access = private)
        
        
        function [buffer_cur_pos] = putMatrixBlockToHessianStorage(this, A, row, col, buffer_cur_pos)
            
            [m, n] = size(A);
            
            buffer_end = m * n + buffer_cur_pos - 1;
            
            this.SparseHessianMatrix_ijv.i_vec(buffer_cur_pos : buffer_end) = kron (ones(1,n), row:row+m-1);
            
            this.SparseHessianMatrix_ijv.j_vec(buffer_cur_pos : buffer_end) = kron (col:col+n-1, ones(1,m));
            
            this.SparseHessianMatrix_ijv.v_vec(buffer_cur_pos : buffer_end) = full (reshape (A, 1, m * n));
            
            buffer_cur_pos = buffer_end + 1;
            
        end        
        

        
        
        function [flcost, nlcost] = SparseNonlinearLeastSquaresOneIteration (this,  lambda)
            vsum = 0;
            this.Tmpbb = 0;
            for ii = 1 : this.numPointClouds
                vsum = vsum + nnz(this.PointClouds(ii).Vis);
            end
            if this.dim == 2
                rot_dim = 1;
                pose_dim = 3;
                const_dim_ones = [1, 1];
                n_elements = 9 * this.numPointClouds + this.dim*this.numPoints + 2 * vsum* this.dim*pose_dim;
            end
            if this.dim == 3
                rot_dim = 3;
                pose_dim = 6;
                const_dim_ones = [1, 1, 1];
                n_elements = 36 * this.numPointClouds + this.dim*this.numPoints + 2 * vsum* this.dim*pose_dim;
            end
            this.SparseHessianMatrix_ijv.i_vec = zeros(1, n_elements);
            this.SparseHessianMatrix_ijv.j_vec = zeros(1, n_elements);
            this.SparseHessianMatrix_ijv.v_vec = zeros(1, n_elements);
            this.SparseRhsVector = zeros(pose_dim*this.numPointClouds + this.dim*this.numPoints,  1);
            
            buffer_cur_pos = 1;
            
            col_offset = this.numPointClouds * pose_dim;
            
            vis_vec_sum_dim = zeros(1, this.dim*this.numPoints);
            vec_rhs_sum = zeros(this.dim*this.numPoints, 1);
            
  
            for kk = 1 : this.numPointClouds                                
                R = this.mPoses{kk}(:, 1:this.dim);
                t = this.mPoses{kk}(:, this.dim+1);
                row_id = pose_dim*kk-pose_dim+1;
                RD = R * full(this.PointClouds(kk).Data(:, logical(this.PointClouds(kk).Vis)));                
                BlockABAB = zeros(pose_dim, pose_dim);
                
                DiffMatrix = (RD + t) - this.mShape(:, logical(this.PointClouds(kk).Vis));
                hk = zeros(pose_dim ,1);
                
                for ii = 1 : size(RD, 2)
                    if this.dim == 2
                        J = [0 , -1; 1, 0] * RD(:, ii);
                    end
                    if this.dim == 3
                        J = -SO3.Hat(RD(:, ii));
                    end
                    PoseJaocbian = [J, eye(this.dim)];
                    BlockABAB = BlockABAB + PoseJaocbian' * PoseJaocbian;
                    col_id = col_offset + (this.PointClouds(kk).PtID(ii)-1)*this.dim+1;
                    [buffer_cur_pos] = this.putMatrixBlockToHessianStorage(PoseJaocbian', row_id, col_id, buffer_cur_pos);
                    [buffer_cur_pos] = this.putMatrixBlockToHessianStorage(PoseJaocbian, col_id, row_id, buffer_cur_pos);
                    hk = hk + PoseJaocbian' * DiffMatrix(:, ii);
                end
                if (lambda > 0)
                    BlockABAB = lambda * eye(pose_dim);
                end
                [buffer_cur_pos] = this.putMatrixBlockToHessianStorage(BlockABAB, row_id, row_id, buffer_cur_pos);
                vis_vec_sum_dim = vis_vec_sum_dim + kron(this.PointClouds(kk).Vis,  const_dim_ones);  %%% this accumulate to CC diagonal
                %%% this accumulate to rhs
                this.SparseRhsVector(row_id : pose_dim*kk) = hk;
                vec_rhs_sum = vec_rhs_sum  + reshape((R * this.PointClouds(kk).Data + t - this.mShape) .* this.PointClouds(kk).Vis, this.dim*this.numPoints ,1);                
                % for LM
                tmp = norm(DiffMatrix, 'fro');
                this.Tmpbb = this.Tmpbb + tmp*tmp;
            end
            if (lambda > 0)
                vis_vec_sum_dim = vis_vec_sum_dim + lambda;
            end
            this.SparseHessianMatrix_ijv.i_vec(buffer_cur_pos : end) = ((col_offset+1) : (col_offset+length(vis_vec_sum_dim)));
            this.SparseHessianMatrix_ijv.j_vec(buffer_cur_pos : end) =  ((col_offset+1) : (col_offset+length(vis_vec_sum_dim)));
            this.SparseHessianMatrix_ijv.v_vec(buffer_cur_pos : end) = vis_vec_sum_dim;
            this.SparseRhsVector((col_offset+1):end) = vec_rhs_sum;
            
            this.SparseRhsVector = - this.SparseRhsVector;
            
            Hessian = sparse(this.SparseHessianMatrix_ijv.i_vec, this.SparseHessianMatrix_ijv.j_vec, this.SparseHessianMatrix_ijv.v_vec);
            
%             size(Hessian)
%             size(this.SparseRhsVector)

            if (lambda > 0)
                this.UpdateX =  Hessian \ this.SparseRhsVector;
            else
                this.UpdateX =  Hessian((pose_dim+1):end, (pose_dim+1):end) \ this.SparseRhsVector((pose_dim+1):end);                
            end
            
            % update poses
            if (lambda == 0)
               init_index = 2; 
            else
               init_index = 1;
            end            
            idx_cnt = 1;           
            for ii = init_index : this.numPointClouds
                deltaX = this.UpdateX( idx_cnt : (idx_cnt+pose_dim-1) );
                idx_cnt = idx_cnt + pose_dim;
                deltaX_rota = deltaX(1:rot_dim);
                if this.dim == 2
                    deltaR = [cos(deltaX_rota),  -sin(deltaX_rota);  sin(deltaX_rota),  cos(deltaX_rota)];
                end
                if this.dim == 3
                    deltaR = SO3.Exp(deltaX_rota);
                end
                R = this.mPoses{ii}(:, 1:this.dim);
                t = this.mPoses{ii}(:, this.dim+1);
                this.mPoses{ii} = [deltaR*R,  t+deltaX((rot_dim+1) : end)];
            end
            % update reference shape;
            this.mShape =  this.mShape - reshape(this.UpdateX(idx_cnt:end),  this.dim, this.numPoints);

            
            if (lambda > 0)
                xAAx = this.UpdateX' * (Hessian - lambda * speye(length(Hessian))) * this.UpdateX;
                xAb = this.UpdateX' * this.SparseRhsVector;
            end
            if (lambda == 0)
                xAAx = this.UpdateX' * Hessian((pose_dim+1):end, (pose_dim+1):end) * this.UpdateX;
                xAb = this.UpdateX' * this.SparseRhsVector((pose_dim+1):end);
            end
            bb = this.Tmpbb;
            flcost = xAAx - 2*xAb + bb;
            nlcost = sqrt(bb/vsum);
            flcost = sqrt(flcost/vsum);

        end


        function [transform_points_cell, transform_vis_cell, refShape] = runTrainingTestData(this, trainingPointsFlag, testingPointsFlag)
            AllPointClouds = this.PointClouds;
            AllNumPoints = this.numPoints;
            this.PointClouds = struct ( 'Data', {}, 'Vis', {},  'PtID', {});
            this.numPointClouds = 0;
            for ii = 1 : length(AllPointClouds)
                data = AllPointClouds(ii).Data(:,  logical(trainingPointsFlag));
                vis = AllPointClouds(ii).Vis(:,  logical(trainingPointsFlag));
                this.addPointCloud (data, vis);
            end
            this.optimize();
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

    
    
    
    
    
    properties (Access = private)
        
        mSqrtLambda = []
        
        EuclideanTransformationSign = []
        
    end
    
    
    
    % methods used to compute Lambda and poses
    methods (Access = private)
        
        
        function initalizePoseByAffineGPA (this)
            Q = zeros(this.numPoints, this.numPoints);
            vis = zeros(1, this.numPoints);
            for ii = 1 : this.numPointClouds
                Pts = [this.PointClouds(ii).Data;  ones(1, this.numPoints) ];
                Pts = Pts .*  this.PointClouds(ii).Vis;
                Q = Q + Pts' * pinv (full(Pts * Pts')) * Pts;
                vis = vis + this.PointClouds(ii).Vis;
                %check_norm = norm(Pts' * pinv (Pts * Pts') * Pts * ones(this.numPoints, 1) - this.PointClouds(ii).Vis')
            end
            Q = diag(vis) - Q;
            if norm(Q * ones(this.numPoints, 1), 'fro') > 1e-6
                fprintf(2, 'Fatal Error. Q*1=0 is violated ||Q*1|| = %.15f \n',   norm(Q * ones(this.numPoints, 1), 'fro'));  %return;
            end
            % Q * 1 = 0.
            Q = Q + (this.numPointClouds/this.numPoints); % soft regularization
            [EV, ~, ~] = this.byEigen (Q, this.dim, 'smallest'); % descending order
            this.mShape =  EV';
            % initialize affine transformations
            for ii = 1 : this.numPointClouds
                Pts = [this.PointClouds(ii).Data .* this.PointClouds(ii).Vis;  ones(1, this.numPoints) ];
                this.mPoses{ii} = full(this.mShape * Pts') * pinv (full(Pts * Pts'));
            end
            % solve for rigid tranformations.  Lambda -> poses
            this.InitializeLambda();
            this.InitializePoses();
            this.mShape = diag(this.mSqrtLambda) * this.mShape;
        end
        
        
        function InitializePoses (this)
            if abs(sum(this.EuclideanTransformationSign)) ~= this.numPointClouds
                fprintf(2, 'Consistent optimal rigid transformation to the reference point-cloud DO NOT exist. \n');
                sign_of_euclidean_transformations = this.EuclideanTransformationSign
            end
            for ii = 1 : this.numPointClouds
                Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));
                % Qt = this.mShape(:, logical(this.PointClouds(ii).Vis));
                Qt = full(this.transformPoints(Pt, ii));
                meanQt = mean(Qt, 2);
                meanPt = mean(Pt, 2);
                R = (Pt - meanPt) * (Qt - meanQt)' * diag(this.mSqrtLambda);
                [U, ~, V] = svd(R); det_VU = det(V*U'); % this has the same sign as det(R). If det(R) > 0, then det(V'U) > 0
                diagS = eye(this.dim);  diagS(end, end) = det_VU;
                R = V * diagS * U';
%                 if (det_VU < 0)
%                     fprintf(2, 'Warning @Det(VU) < 0 in orthogonal Procrustes \n');
%                 end
                t = - (R * meanPt - diag(this.mSqrtLambda) * meanQt);
                this.mPoses{ii} = [R, t];
            end
        end
        

        function Lambda = InitializeLambda (this)
            this.EuclideanTransformationSign = zeros(1, this.numPointClouds);
            AMatrix = zeros(this.dim * this.dim * this.numPointClouds, this.dim);
            bVector = zeros(this.dim * this.dim * this.numPointClouds, 1);
            n_elements = this.dim * this.dim;
            inv_lambda = zeros(this.dim, 1);
            for ii = 1 : this.numPointClouds
                Pt = full(this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis)));
                % Qt = this.mShape(:, logical(this.PointClouds(ii).Vis));
                Qt = full(this.transformPoints(Pt, ii));
                zcQt = Qt - mean(Qt, 2);
                zcPt = Pt - mean(Pt, 2);
                Lt = full(zcQt * zcPt' * pinv(full(zcPt * zcPt')));
                this.EuclideanTransformationSign(ii) = sign(det(Lt));  % the sign of affine transformation
                %  || L * L' - Inv(Lambda) ||_F = || diag( L * L' - Inv(Lambda) ) ||_2
                tmp = Lt * Lt';
                if (length(Lt) == 2)
                    inv_lambda = inv_lambda + [tmp(1,1); tmp(2,2)];
                end
                if (length(Lt) == 3)
                    inv_lambda = inv_lambda + [tmp(1,1); tmp(2,2); tmp(3,3)];
                end
            end
            inv_lambda = inv_lambda / this.numPointClouds;
            Lambda = 1 ./ inv_lambda;
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
%                 Lt = full(zcQt * zcPt' * pinv(full(zcPt * zcPt')));
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
%                 this.mPoses{ii}(1,:) = - this.mPoses{ii}(1,:);
%                 this.EuclideanTransformationSign = - this.EuclideanTransformationSign;
%             end
%         end


        
        
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
    
    
    
    
end

