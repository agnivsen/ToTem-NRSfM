classdef GPA_Functor < handle

    properties (Access = public)
        
        dim = 0
        
        numPoints = 0
        
        mShape = []
        
        mPoses = {}
        
        numPointClouds = 0
        
        rsd_ref = 0
        
        rsd_dat = 0
        
        time = 0
        
        dimFeatureSpace = -1;
        
        verbose = 'off'

    end
    
    properties (Access = private)
        
        method
        
        mDataCell
        
        mVisVecCell
        
        A = {}
        
        a = {}
        
        maxNum = 10000000

    end

    methods( Access = public)
        
        function this = GPA_Functor (method_str)
            
            addpath('DefGPA');
            addpath('SGPAv1.0', 'SGPAv1.0/util');
            
            this.numPointClouds = 0;
            
            this.method = 'EUC-ALL';
            
            if (nargin > 0)
                
                this.method = method_str;
            
            end
            
           % fprintf(1, '%s\n', this.method);
            
        end
           
    end
    
    
    
    methods( Access = public)
        
        function addPointCloud (this, dataMatr, visVec )
            
            cnt = this.numPointClouds + 1;
            
            this.mDataCell{cnt} = dataMatr;
            
            this.mVisVecCell{cnt} = ( sum (abs(dataMatr) < this.maxNum) == size(dataMatr, 1) );
            
            if (nargin == 3)
                
                this.mVisVecCell{cnt} = visVec;
                
            end
            
            this.numPointClouds = cnt;
            
            this.dim = size(dataMatr, 1);
            
            this.numPoints = size(dataMatr, 2);
            
        end
        
        
        function addDataByArray(this, DataCell, VisVecCell)
            for ii = 1 : length(DataCell)
                this.addPointCloud(DataCell{ii}, VisVecCell{ii});
            end            
        end

        
        function run (this, smoothParam)
                        
            options.method = this.method;
            
            options.verbose = this.verbose;
            
            gpaResultDataSpaceCost = gpa(this.mDataCell, this.mVisVecCell, options);
            
            this.time(1) = gpaResultDataSpaceCost.time(1);
            
            % data space cost is invariant to a global gauge transformation
            % transform gauge to the coordinate frame of the first datum shape
            G = gpaResultDataSpaceCost.A{1};
            
            g = gpaResultDataSpaceCost.a{1};
            
            S = G * gpaResultDataSpaceCost.S + g;
            
            meanS = mean(S,2);
            
            this.mShape = S - meanS;
            
            for ii = 1 : this.numPointClouds
                
                this.A{ii} = gpaResultDataSpaceCost.A{ii} * inv(G);
                
                this.a{ii} = - gpaResultDataSpaceCost.A{ii} * inv(G) * (g - meanS) + gpaResultDataSpaceCost.a{ii};
                
            end

            for ii = 1 : this.numPointClouds
                
                invA = inv(this.A{ii});
                
                this.mPoses{ii} = [invA,  - invA*this.a{ii}];
                
            end
            
        end

        
        function [rsd_ref, rsd_dat] = rsdError(this)
            
            rsd_ref = 0;  rsd_dat = 0;  v = 0;
            
            %biggest =0; iddd = 1;
            
            for ii = 1 : this.numPointClouds
                
                v = v + nnz(this.mVisVecCell{ii});
                
                % reference space cost
                
                Diff = (inv(this.A{ii}) * this.mDataCell{ii} - inv(this.A{ii}) * this.a{ii}) - this.mShape;
                
                dr = norm(Diff .* this.mVisVecCell{ii}, 'fro');
                
                rsd_ref = rsd_ref + dr * dr;
                
                % datum space cost
                
                Diff = this.A{ii} * this.mShape + this.a{ii} - this.mDataCell{ii};
                
                dd = norm(Diff .* this.mVisVecCell{ii}, 'fro');
                
                rsd_dat = rsd_dat + dd * dd;
                
%                 if (dd * dd > biggest)
%                     
%                     biggest = dd * dd;
%                     
%                     iddd = ii;
%                     
%                     fprintf(2, 'val = %f,  idx = %d\n', biggest, iddd);
%                     
%                 end
                
            end
            
            rsd_ref = sqrt(rsd_ref/v);
            
            rsd_dat = sqrt(rsd_dat/v);
            
            this.rsd_ref = rsd_ref;  this.rsd_dat = rsd_dat;
            
        end
        
        
        function pImageArray = transformPoints (this, pValArray, cloudIndex )

            pImageArray = inv(this.A{cloudIndex}) * pValArray - inv(this.A{cloudIndex}) * this.a{cloudIndex};
                        
        end
        
        
        function ptCloud = generatePointCloud (this, cloudIndex)
            
            ptCloud = this.A{cloudIndex} * this.mShape + this.a{cloudIndex};
            
        end
        
    end

end