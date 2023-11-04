classdef WarpProcrustesAnalysis < handle
    
    
    % input data
    properties (Access = public)
        
        PointClouds = struct ( 'Data', {}, 'Vis', {} )
        
        numPointClouds = 0
        
        numPoints = 0
        
        warpModelFlag = 'TPS'
        
        smoothParam = 1 % param for regularizations
        
        translationParam = 1 % param for zero translations
        
        dimFeatureSpace = 10
        
        maxNum = 10000  % the celing number, should never be reached.
        
        fileExt = 'txt'
        
        useZeroCenteredData = 1
        
    end
    
    
    % result
    properties (Access = public)
        
        mShape = []
        
        mWeight = []
        
    end
    
    
    % internal feature space data
    properties (Access = private)
       
        FeatureClouds = struct ( 'FeatureMatr', {}, 'VisMatr', {}, 'Regularization', {} )
        
        numVisPoints = []
        
    end
    
    
    % sparse matrices
    properties (Access = private)
        
        VisVecGlob = []
        
        Dt = struct ('i', [], 'j', [], 'v', [])
        
        Nt = struct ('i', [], 'j', [], 'v', [])

        Z = struct ('i', [], 'j', [], 'v', [])
        
        Ft = struct ('i', [], 'j', [], 'v', [])

        InvXi = struct ('i', [], 'j', [], 'v', [])
        
    end
    
    
    % warp models
    properties (Access = private)
        
       TPS = struct ( 'IntrisicMatr', {}, 'ControlPoints', {}, 'InternalSmoothParam', {} )
        
    end
    
    
    

    
    
    methods (Access = public)
 
        % class constructor
        function this = WarpProcrustesAnalysis ( ceilingNum )
            
            if (nargin > 0)
            
                this.maxNum = ceilingNum;
            
            end
            
        end
        
        
        function addPointCloud (this, dataMatr )
           
            cnt = this.numPointClouds + 1;

            this.PointClouds(cnt).Data = dataMatr;
            
            this.PointClouds(cnt).Vis = ( sum (abs(dataMatr) < this.maxNum) == 3 );
            
            this.numPointClouds = cnt;
            
        end
        
        
        function run (this)
            
            this.rePositionData();
            
            tic
            
            this.initWarpModel ();
            
            usedTime = toc;
            disp(['time used: initialize warp model:  ', num2str(usedTime)]);
            tic
            
            this.toFeatureCloud();

            usedTime = toc;
            disp(['time used: mapping to feature cloud:  ', num2str(usedTime)]);
            tic
            
            this.setDt();
            
            this.setNt();

            this.setInvXi();
            
            usedTime = toc;
            disp(['time used: compute vectorized sparse matrices:  ', num2str(usedTime)]);
            tic
            
            mDt = sparse (this.Dt.i, this.Dt.j, this.Dt.v);
            
            mNt = sparse (this.Nt.i, this.Nt.j, this.Nt.v);
            
            mInvXi = sparse (this.InvXi.i, this.InvXi.j, this.InvXi.v);

            usedTime = toc;
            disp(['time used: construct sparse matrices:  ', num2str(usedTime)]);
            tic            

            PMatr = mNt * mNt' - mNt * mDt' * mInvXi * mDt * mNt';
            
            usedTime = toc;
            disp(['time used: construct matrix for decomposition:  ', num2str(usedTime)]);
            tic      
            
            [EV, ~, EA] = this.byEigen (full(PMatr), 3);
            
            usedTime = toc;
            disp(['time used: solve the shape by Eigen Value Decomposition:  ', num2str(usedTime)]);
            tic
            
            eigVal = this.getEigenValue();  [sVal, sIndex] = sort(eigVal, 'ascend');
            
            this.mShape(sIndex, :) =  diag(sqrt(sVal)) * EV';
            
            this.mWeight = mInvXi * mDt * mNt' * this.mShape';
            
            usedTime = toc;
            disp(['time used: construct the shape & weight:  ', num2str(usedTime)]);
            
            
            this.plotCloud3 (this.PointClouds(1).Data, EV', this.mShape, 'Original data shape noised', 'Unscaled Reference', 'Scaled Reference');
            
            Eigenvalue_FirstN_LastN = [EA(1:10), EA(end-9:end)]
            
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
            
            this.numPoints = size(dataMatr, 2);
            
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
        
        
        function pImage = WarpPoint ( pVal, cloudIndex )
        
            fImage = this.fMapping ( pVal, cloudIndex );
            
            fvec = ( (cloudIndex-1) * this.dimFeatureSpace + 1 ) : ( cloudIndex * this.dimFeatureSpace ); 
            
            pImage = this.mWeight(fvec, :)' *  fImage;
            
        end
        
    end
    
    
    
    % methods used to constructed internal matrices
    methods (Access = private)
  
        
        function toFeatureCloud (this)
            
           for ii = 1 : this.numPointClouds
              
               % construct data mapping matrix
               
               this.FeatureClouds(ii).FeatureMatr = zeros ( this.dimFeatureSpace, this.numPoints );
                              
               for jj = 1 : this.numPoints
                 
                   this.FeatureClouds(ii).FeatureMatr(:, jj) = this.fMapping ( this.PointClouds(ii).Data(:, jj), ii );
               
               end
               
               % construct visiblity matrix
               
               this.numVisPoints(ii) = nnz (this.PointClouds(ii).Vis);
               
               T = sparse (1 : this.numPoints, 1 : this.numPoints, this.PointClouds(ii).Vis);
               
               this.FeatureClouds(ii).VisMatr = T(:, any(T));
               
               % construct bending energy matrix
               
               this.FeatureClouds(ii).Regularization = this.getRegularization (ii);
               
           end
            
        end
        
        
        
        function setDt (this)
            
            eleNonzeros = this.dimFeatureSpace * sum (this.numVisPoints);
            
            this.Dt.i =  ones (1, eleNonzeros);
            
            this.Dt.j =  ones (1, eleNonzeros);
            
            this.Dt.v = zeros (1, eleNonzeros);
            
            pt = 1; cnt_i = 1; cnt_j = 1;
            
            for ii = 1 : this.numPointClouds
                
                eleNum = this.dimFeatureSpace * this.numVisPoints(ii);
                
                [si, sj, sv] = this.putMatrixBlock (this.FeatureClouds(ii).FeatureMatr * this.FeatureClouds(ii).VisMatr, cnt_i, cnt_j);
                
                this.Dt.i(pt:pt+eleNum-1) = si;
                
                this.Dt.j(pt:pt+eleNum-1) = sj;
                
                this.Dt.v(pt:pt+eleNum-1) = sv;
                
                pt = pt + eleNum;
                
                cnt_i = cnt_i + this.dimFeatureSpace;
                
                cnt_j = cnt_j + this.numVisPoints(ii);
                
            end

        end
        
        
        
        function setNt (this)
            
            eleNonzeros = this.numPoints * sum (this.numVisPoints);
            
            this.Nt.i =  ones (1, eleNonzeros);
            
            this.Nt.j =  ones (1, eleNonzeros);
            
            this.Nt.v = zeros (1, eleNonzeros);
            
            pt = 1; cnt_i = 1; cnt_j = 1;
            
            for ii = 1 : this.numPointClouds
                
                eleNum = this.numPoints * this.numVisPoints(ii);
                
                [si, sj, sv] = this.putMatrixBlock (this.FeatureClouds(ii).VisMatr, cnt_i, cnt_j);
                
                this.Nt.i(pt:pt+eleNum-1) = si;
                
                this.Nt.j(pt:pt+eleNum-1) = sj;
                
                this.Nt.v(pt:pt+eleNum-1) = sv;
                
                pt = pt + eleNum;
                
                cnt_j = cnt_j + this.numVisPoints(ii);
                
            end            
            
        end
  
        
        
        function setInvXi (this)
            
            eleNonzeros = this.dimFeatureSpace * this.dimFeatureSpace * this.numPointClouds;
            
            this.InvXi.i =  ones (1, eleNonzeros);
            
            this.InvXi.j =  ones (1, eleNonzeros);
            
            this.InvXi.v = zeros (1, eleNonzeros);
            
            pt = 1; cnt_i = 1; cnt_j = 1;
            
            for ii = 1 : this.numPointClouds
                
                eleNum = this.dimFeatureSpace * this.dimFeatureSpace;
                
                DtMatr = this.FeatureClouds(ii).FeatureMatr * this.FeatureClouds(ii).VisMatr;
                
                ZMatr = this.FeatureClouds(ii).Regularization;
                
                FtMatr = this.FeatureClouds(ii).FeatureMatr * this.VisVecGlob;
                
                tmp = DtMatr * DtMatr' + this.smoothParam * ZMatr' * ZMatr + this.translationParam * FtMatr * FtMatr';
               
                tmpMatr = inv( tmp );
                
                [si, sj, sv] = this.putMatrixBlock (tmpMatr, cnt_i, cnt_j);
                
                this.InvXi.i(pt:pt+eleNum-1) = si;
                
                this.InvXi.j(pt:pt+eleNum-1) = sj;
                
                this.InvXi.v(pt:pt+eleNum-1) = sv;
                
                pt = pt + eleNum;
                
                cnt_i = cnt_i + this.dimFeatureSpace;
                
                cnt_j = cnt_j + this.dimFeatureSpace;
                
            end
            
        end
        
        
        
        function setVisVecGlob (this)
            
            tmp = ones (1, this.numPoints);
            
            for ii = 1 : this.numPointClouds
                
                tmp = tmp .* this.PointClouds(ii).Vis;
                
            end
            
            this.VisVecGlob = tmp';
            
        end
        
        
        function rePositionData (this)
            
            this.setVisVecGlob();
            
            % Center the data to zero
            if this.useZeroCenteredData
                
                for ii = 1 : this.numPointClouds
                    
                    gAnchor = (this.PointClouds(ii).Data * this.VisVecGlob)/nnz(this.VisVecGlob);
                    
                    this.PointClouds(ii).Data = this.PointClouds(ii).Data - gAnchor;
                    
                end
                
            end
            
        end
        
        
        function v = getEigenValue (this)
            
            eigVal = zeros (3, this.numPointClouds);
            
            for ii = 1 : this.numPointClouds
                
                tmp = this.PointClouds(ii).Data * this.PointClouds(ii).Data';
                
                eigVal(:, ii) = sort ( eig ( (tmp+tmp')/2 ) );
                
            end
            
            v = sum(eigVal, 2) / this.numPointClouds;
            
        end
        
        
    end
    
    
  
    

    
    
    % utility functions
    methods (Access = public, Static = true)
        
        % a function used to construct a sparse matrix by given sublocks
        % this function put a matrix A at the position (row, col) of a sparse matrix.
        % use " sparse (i, j, v) " to check the resulting sparse matrix
        function [i, j, v] = putMatrixBlock (A, row, col)
            
            [m, n] = size(A);
            
            v = full (reshape (A, 1, m * n));
            
            i = kron (ones(1,n), row:row+m-1);
            
            j = kron (col:col+n-1, ones(1,m));
            
        end
        
        
        % Thin-Plate-Spline basis function for 3D case
        % it's a radial distance function
        % dp = the differnece of two points
        function dist = TPS_BasisFunc ( dp )
            
            r = norm (dp);
            
            dist = - r;
            
        end
        
        
        function plotCloud3 (Cloud1, Cloud2, Cloud3, figTitle1, figTitle2, figTitle3)
            
            fig = figure;
            fig.Position = [50, 800, 1800, 800];
            
            subplot(1,3,1)
            
            plot3(Cloud1(1,:), Cloud1(2,:), Cloud1(3,:), 'r.');
            axis equal;
            title (figTitle1);
            
            subplot(1,3,2)
            plot3(Cloud2(1,:), Cloud2(2,:), Cloud2(3,:), 'b.');
            axis equal;
            title (figTitle2);
            
            subplot(1,3,3)
            plot3(Cloud3(1,:), Cloud3(2,:), Cloud3(3,:), 'g.');
            axis equal;
            title (figTitle3);
            
            sgtitle('Shape Estimations: All data shapes are generated by noising the original data shape, then adding rotations.');
            
        end
        
    end
    
    
    
    % methods to solve optimziation problem: 
    %  min || K X ||_F
    %  s.t.  X^T X = I
    methods (Access = private, Static = true)
        
        % K, the matrix to decompose,
        % d, the number of eigen vectors required.
        function [V, eigValEigen, eigenValueAll] = byEigen (K, d)
            
            offset = 0;

            [Q, D ] = eig( (K+K')/2 );
            
            eigVal = diag(D);
            
            [eigVal, Index] = sort(eigVal, 'descend');  Q = Q (:, Index);

            V = Q (:, end-d+1-offset:end-offset);
            
            eigenValueAll = eigVal;
            
            eigValEigen =  eigVal (end-d+1-offset:end-offset);
            
            CheckEigenVectorOrthonormal = V' * V

        end

        
        function [V, eigValSVD, singularValueAll] = bySVD (K, d)
            
            offset = 0;

            [~, D, V] = svd(K);
            
            V = V(:, end-d+1-offset:end-offset);
            
            singularValueAll = diag(D);
            
            eigValSVD = singularValueAll (end-d+1-offset:end-offset).^2;

        end
 
        
    end
    
    
    
    % These two methods are differnt for different warp models.
    % implement with flags as an option
    methods (Access = private)
        
        function initWarpModel (this)
            
            if strcmp (this.warpModelFlag, 'TPS')
                
                for ii = 1 : this.numPointClouds
                    
                    % set control points of TPS

                    min_xyz = min (this.PointClouds(ii).Data, [], 2);
                    max_xyz = max (this.PointClouds(ii).Data, [], 2);
                    spn_xyz = max_xyz - min_xyz;
                    
                    this.TPS(ii).ControlPoints = [ min_xyz(1) + rand(1, this.dimFeatureSpace) * spn_xyz(1);
                                                   min_xyz(2) + rand(1, this.dimFeatureSpace) * spn_xyz(2);
                                                   min_xyz(3) + rand(1, this.dimFeatureSpace) * spn_xyz(3); ];
                    
                    % set internal smoothing parameter of TPS
                    
                    this.TPS(ii).InternalSmoothParam = 0;
                    
                    % compute the intrisic matrix of TPS
                    
                    K = zeros (this.dimFeatureSpace, this.dimFeatureSpace);
                    
                    for it_r = 1 : this.dimFeatureSpace
                        
                        for it_c = it_r+1 : this.dimFeatureSpace
                            
                            K (it_r, it_c) = this.TPS_BasisFunc( this.TPS(ii).ControlPoints(:, it_r) - this.TPS(ii).ControlPoints(:, it_c) );
                            
                        end
                        
                    end
                    
                    K = K + K';
                    
                    for it_d = 1 : this.dimFeatureSpace
                        
                        K (it_d, it_d) = this.TPS(ii).InternalSmoothParam;
                        
                    end
                    
                    InvK = inv (K);
                    
                    C = [ this.TPS(ii).ControlPoints; ones(1, this.dimFeatureSpace) ];
                    
                    H = inv(C * InvK * C') * C * InvK;
                    
                    this.TPS(ii).IntrisicMatr = [ InvK - InvK * C' * H;  H ];
                    
                end
                
                return
                
            end
            
        end
        
        
        
        % mapping R^3 -> R^l : point to feature space
        function fImage = fMapping (this, pVal, cloudIndex)
           
            if strcmp (this.warpModelFlag, 'TPS')
                
                nlv = ones (this.dimFeatureSpace + 3 + 1, 1);
                
                dmatr = this.TPS(cloudIndex).ControlPoints - pVal * ones (1, this.dimFeatureSpace);
                
                for k = 1 : this.dimFeatureSpace
                    
                    nlv (k) = this.TPS_BasisFunc (dmatr(:, k));
                    
                end
                
                nlv (this.dimFeatureSpace+1 : this.dimFeatureSpace+3) = pVal;
                
                fImage = this.TPS(cloudIndex).IntrisicMatr' * nlv;
                
                return
                
            end
            
        end
        
        
        % set the bending energy matrix for each FeatureCloud
        function Z = getRegularization (this, cloudIndex)
            
            if strcmp (this.warpModelFlag, 'TPS')
                
                M = 8 * pi * this.TPS(cloudIndex).IntrisicMatr (1:this.dimFeatureSpace, :);
                
                [V, D]= eig ( (M+M')/2 );  [Q, ~] = qr (V);
                
                for ii = 1 : size(D,1)
                    if  abs (D(ii,ii)) < 1e-12
                        D(ii, ii) = 0;
                    else
                        D(ii, ii) = sqrt (abs (D(ii,ii)));
                    end
                end
                
                Z = Q * D * Q';
                
                if   norm  (Z'*Z - M, 'f') > 1e-12
                    
                    fprintf(2, 'Bending energy matrix decomposition fails.');
                
                end
                
                return
                
            end
            
        end
                
    end
    
    
end

