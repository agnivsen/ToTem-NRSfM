classdef DGPA_Affine < handle
    
    
    % input data
    properties (Access = public)
        
        dim = 0
        
        PointClouds = struct ( 'Data', {}, 'Vis', {}, 'TransPrior', {} )
        
        numPointClouds = 0
        
        numPoints = 0
        
        maxNum = 1000000  % the celing number, should never be reached.
        
        fileExt = 'txt'
        
        useZeroCenteredData = true
        
    end
    
    
    % result
    properties (Access = public)
        
        mShape = []
        
        AffineTransform = struct ('Affine', {}, 'Trans', {})
        
        rsd_ref = 0
        
        rsd_dat = 0
        
        rigidity = 0
        
        time = []

    end
    
    
    % internal matrices
    properties (Access = private)
        
        VisVecGlob = []
        
        P = []
        
    end
    
    


    
    methods (Access = public)
 
        % class constructor
        function this = DGPA_Affine ( ceilingNum )
            
            if (nargin > 0)
            
                this.maxNum = ceilingNum;
            
            end
            
        end
        
        
        function addPointCloud (this, dataMatr, visVec )

            cnt = this.numPointClouds + 1;

            this.PointClouds(cnt).Data = dataMatr;
            
            this.PointClouds(cnt).Vis = ( sum (abs(dataMatr) < this.maxNum) == size(dataMatr, 1) );
            
            if (nargin == 3)
                
                this.PointClouds(cnt).Vis = visVec;
                
            end
            
            this.numPointClouds = cnt;
            
            this.dim = size(dataMatr, 1);
            
            this.numPoints = size(dataMatr, 2);
            
        end
        
        
        function run (this)
            
            tic;
            
            this.rePositionData();
            
            this.setP();
            
            [EV, eigP, EA] = this.byEigen (this.P, this.dim, 'smallest'); % descending order
            
            eigVal = this.getSingularValue(); % ascending order
            
            this.mShape =  diag(sqrt(eigVal)) * EV';
            
            % check valid rotation, and eliminate reflection
            sR = this.ScaledOrhoProcurestes (this.PointClouds(1).Data, this.mShape);
            if (det(sR) < 0)
                this.mShape(1,:) = - this.mShape(1,:);
            end
            
            this.setAffineTransform();
                   
            this.time(1) = toc;
            
            %disp(['Cost function: ', num2str(eigP' * eigVal)]);
            
%             Eigenvalue_FirstN_LastN = [EA(1:10), EA(end-9:end)];
         
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
        

         
        function pImageArray = transformPoints (this, pValArray, cloudIndex )
            
            pImageArray = pValArray + this.PointClouds(cloudIndex).TransPrior;
            
            pImageArray = this.AffineTransform(cloudIndex).Affine * pImageArray;
            
        end
        
        
        
        function [rsd_ref, rsd_dat] = rsdError (this)
            
            rsd_ref = 0; rsd_dat = 0; v = 0;
            
            for ii = 1 : this.numPointClouds
                
                v = v + nnz(this.PointClouds(ii).Vis);
                
                Diff_ref = this.mShape - this.AffineTransform(ii).Affine * this.PointClouds(ii).Data;
                
                Diff_dat = inv(this.AffineTransform(ii).Affine) * this.mShape -  this.PointClouds(ii).Data;
                
                dr = norm( Diff_ref .* this.PointClouds(ii).Vis, 'fro');
                
                dd = norm( Diff_dat .* this.PointClouds(ii).Vis, 'fro');
                
                rsd_ref = rsd_ref + dr * dr;
                
                rsd_dat = rsd_dat + dd * dd;
                
            end
            
            rsd_ref = sqrt(rsd_ref/v);
            
            rsd_dat = sqrt(rsd_dat/v);
            
            this.rsd_ref = rsd_ref;  this.rsd_dat = rsd_dat;
            
        end

        
        
        function hfig = plotPointCloud (this)
            
            n = this.numPointClouds;
        
            hfig = figure;
            hfig.Position = [50, 800, 1800, 800];
          
            DT = delaunayTriangulation(this.mShape');
            
            for ii = 1 : n
                
                subplot(2, ceil(n/2),ii);
                
                if this.dim == 3
                
                    [TRI, ~] = convexHull (DT);
                    
                    trisurf(TRI, this.PointClouds(ii).Data(1,:), this.PointClouds(ii).Data(2,:), this.PointClouds(ii).Data(3,:), 'LineWidth', 1.5, 'EdgeColor', 'r', 'FaceColor', 'b');
                
                end
                
                if this.dim == 2
                    
                    TRI = DT.ConnectivityList;
                    
                    triplot(TRI, this.PointClouds(ii).Data(1,:), this.PointClouds(ii).Data(2,:), '-b.', 'MarkerSize', 15);
                    
                end
                
                axis equal;
                
                title (num2str(ii));
            
            end
               
            sgtitle('Datum Shapes');
            
         end
        
        
         function hfig = plotRefShape (this)
             
             hfig = figure;
             title('Optimized Shapes');
             hold on
             
             DT = delaunayTriangulation(this.mShape');
             
             for ii = 1 : this.numPointClouds
                 
                 D = this.AffineTransform(ii).Affine * this.PointClouds(ii).Data;
                 
                 if this.dim == 3
                     
                     [TRI, ~] = convexHull (DT);
                     
                     trisurf (TRI, D(1,:), D(2,:), D(3,:), 'LineWidth', 1.0, 'EdgeColor', 'r', 'FaceColor', 'b');
                     
                 end
                 
                 if this.dim == 2
                     
                     TRI = DT.ConnectivityList;
                     
                     triplot (TRI, D(1,:), D(2,:),  '-b.',  'MarkerSize', 15);
                     
                 end
                 
             end
             
             if this.dim == 3
                 
                 [TRI, ~] = convexHull (DT);
                 
                 trisurf (TRI, this.mShape(1,:), this.mShape(2,:), this.mShape(3,:), 'LineWidth', 1.5, 'EdgeColor', 'r', 'FaceColor', 'b');
                 
                 view(3);
                 
             end
             
             if this.dim == 2
                 
                 TRI = DT.ConnectivityList;
                 
                 triplot (TRI, this.mShape(1,:), this.mShape(2,:),  '-g.', 'LineWidth', 1.5,  'MarkerSize', 20);
                 
             end
             
             axis equal;
             hold off
             
         end

         
         
         function rig = rigidityScore (this)
             
             rig = zeros(this.numPointClouds, 1);
                          
             for ii = 1 : this.numPointClouds
                 
                 D = this.PointClouds(ii).Data(:, logical(this.PointClouds(ii).Vis));
                 
                 Dtrans = this.AffineTransform(ii).Affine * D;
                 
                 rig (ii) = this.DistScore (D, Dtrans);
                 
             end
             
             this.rigidity = 1 - rig;
             
         end
         

    end
    
    
    
    % methods used to constructed internal matrices
    methods (Access = private)
  
        

        function setP (this)
            
            this.P = zeros(this.numPoints, this.numPoints);
            
            for ii = 1 : this.numPointClouds
               
                this.P = this.P + this.PointClouds(ii).Data' * inv (this.PointClouds(ii).Data * this.PointClouds(ii).Data') * this.PointClouds(ii).Data;
       
            end
            
            vis = zeros(1, this.numPoints);
            
            for ii = 1 : this.numPointClouds
                
                vis = vis + this.PointClouds(ii).Vis;
                
            end
            
            this.P = diag(vis) - this.P;
            
        end
        
        
        function setAffineTransform (this)
            
            for ii = 1 : this.numPointClouds
                
                At = inv (this.PointClouds(ii).Data * this.PointClouds(ii).Data') * this.PointClouds(ii).Data * this.mShape';
                
                this.AffineTransform(ii).Affine =  At';
                
                this.AffineTransform(ii).Trans = At' * this.PointClouds(ii).TransPrior;
                
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
                    
                    this.PointClouds(ii).Data = this.PointClouds(ii).Data .* this.PointClouds(ii).Vis;
                    
                    this.PointClouds(ii).TransPrior = - gAnchor;
                    
                end
                
            end
            
        end
        
        
        
        function v = getSingularValue (this)
            
            arrayD = cell(1, this.numPointClouds);
            
            sumVis = zeros(1, this.numPoints);
            
            for ii = 1 : this.numPointClouds
            
                arrayD{ii} = this.PointClouds(ii).Data (:, logical(this.VisVecGlob));
                
                sumVis = sumVis + this.PointClouds(ii).Vis;
            
            end
            
            eigVal = zeros (this.dim, this.numPointClouds);
            
            for ii = 1 : this.numPointClouds
                
                sumD = zeros(this.dim, this.numPoints);
                
                for jj = 1 : this.numPointClouds
                    
                    sRij = this.ScaledOrhoProcurestes (arrayD{jj}, arrayD{ii});
                    
                    sumD = sumD + sRij * this.PointClouds(jj).Data;
                    
                end
                
                vcnt = sumVis - this.PointClouds(ii).Vis;
                
                vcnt(vcnt > 0) = 1 ./ vcnt(vcnt>0);
                
                vcnt = vcnt .* (1 - this.PointClouds(ii).Vis);
                
                D = this.PointClouds(ii).Data .* this.PointClouds(ii).Vis + (sumD - this.PointClouds(ii).Data) .*  vcnt;
                
                tmp = D * D';
                
                eigVal(:, ii) = sort ( eig ( (tmp+tmp')/2 ), 'ascend' );
                
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
        
        
    end
    
    
    
    % methods to solve optimziation problem: 
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
            
            if  norm( V' * V - eye(d) ) > 1e-12
                
                fprintf(2, 'Eigen vectors are not orthnormal!.');
                
            end

        end
        
        
        
        % The solution to the scaled orthogonal Procrustes problem
        % sR = minimize || sR * D1  - D2 ||
        function sR = ScaledOrhoProcurestes (D1, D2)
            
            meanD1 = mean(D1, 2);
            
            meanD2 = mean(D2, 2);
            
            M = (D1 - meanD1) * (D2 - meanD2)';
            
            [U, ~, V] = svd(M);
            
            R = V * U';
            
            if det(R) < 0
                
                fprintf(2, '\nThere exists reflection between solutions! \n');
                
            end
            
            s = trace(R * M) / trace(M);
            
            sR = s * R;
            
        end
 
        
        
        % distance score between two shapes
        % used to evaluate rigidity of transformations
        function distScore = DistScore (D, Dtrans)
            
            DT = delaunayTriangulation (D');
            
            edges = DT.edges();
            
            diff = zeros (size(edges,1), 1);
            
            for jj = 1 : size(edges, 1)
                
                e = edges(jj, :);
                
                distD = norm( D(:, e(1)) - D(:, e(2)) );
                
                distDtrans = norm ( Dtrans(:, e(1)) - Dtrans(:, e(2)) );
                
                diff(jj) =  abs ((distD - distDtrans) / distD);
                
            end
            
            distScore = mean (diff);
            
        end
        
        
    end
    
    
    

    
end

