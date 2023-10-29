classdef PlotSurfaceReconstruction < handle
    properties(Access = private)
        densePointCloud = []; % to be used for surface reconstruction, dimension [N x 3]
        
        sparsePointCloud = []; % optional pointcloud that could be overlaid on reconstructed surface, dimension [n x 3]

        sparsePointSize = 3;
        
        % OPTIONAL variables:
        sparseNormals = []; % sparse pointcloud could also be associated with normals, dimension [n x 3], should have one-to-one correspondence with 'sparsePointCloud'
        normalLength = 0.1; % length of normals while plotting
        
        lightModel = [];
        material = 0; % default material is DULL (0), could be modified to: 1) 'metal' (1) or 2) 'shiny' (2)
        
        nDensePcl = 0; % no. of pt.s in dense pointcloud
        nSparsePcl = 0; % no. of pt.s in sparse pointcloud
        
        hasDensePointcloud = false;
        hasSparsePointcloud = false;
        hasSparseNormals = false;
        
        showAxis = false; % instruct the code to display axes
        
        numLightModel = 0; % number of lights added to the scene
        
        plotPositionAndSize = [10 10 1000 1000];
        
        outputPath = '';
        hasOutputPath = false;

        img = [];
        mask = [];
        extrinsics = [];
        intrinsics = [];
        hasTexture = false;
        hasMask = false;
        textureDensity = 500;% implying 500 x 500 pixels
    end
    
    methods(Access = public)
        function [obj] = PlotSurfaceReconstruction(pointcloud)
            
            assert(size(pointcloud,1) > 3,'As a trivial sanity check, we expect the dense pointcloud to have at least three points, i.e., the minimum requirement for a surface');
            assert(size(pointcloud,2) == 3,'The dense pointcloud needs to be a [N x 3] array');
            
            obj.densePointCloud = pointcloud;
            obj.hasDensePointcloud = true;
        end
        
        function [obj] = addSparsePointcloud(obj, sparsePointcloud, sparseNormals)
            assert(size(sparsePointcloud,2) == 3,'The sparse pointcloud needs to be a [n x 3] array');
            
            obj.sparsePointCloud = sparsePointcloud;
            obj.hasSparsePointcloud = true;
            
            if exist('sparseNormals', 'var')
                assert(size(sparseNormals,2) == 3,'The sparse normals needs to be a [n x 3] array, this is [OPTIONAL]');
                assert(size(sparseNormals,1) == size(sparsePointcloud,1),'The sparse normals needs to be a [n x 3] array, the same dimension as sparse pointcloud. This variable is [OPTIONAL]');
                
                obj.sparseNormals = sparseNormals;
                obj.hasSparseNormals = true;
            end    
        end
        
        function [obj] = setNormalLength(obj, normlLength)
            assert(normlLength > 0, 'The length of normals need to be greater than zero; if you wish to disable normal visualization, do not set sparse normals in <addSparsePointcloud>');
            obj.normalLength = normlLength;
        end
        
        function [obj] = addLightToScene(obj, position, style)
            if exist('style', 'var')
                assert(strcmp(style, 'local') || strcmp(style, 'infinite'),'Only styles available for light are <local> or <infinite> [default]');
            else
                style = 'infinite';
            end
            
            assert(numel(position)==3, 'Position of light needs to be a [1 x 3] array');
            assert(size(position,1)==1, 'Position of light needs to be a [1 x 3] array');
            
            lightElement.style = style;
            lightElement.position = position;
            
            obj.numLightModel = obj.numLightModel + 1;
            
            obj.lightModel{obj.numLightModel} = lightElement;            
        end
        
        function [obj] = setMaterial(obj, materialType)
            assert(materialType==0 || materialType==1 || materialType==2 , 'Acceptible material types are: 0 - dull, 1 - metal, 2 - shiny');
            obj.material = materialType;
        end

        function [obj] = setSparsePointSize(obj, pointSize)
            obj.sparsePointSize = pointSize;
        end

        function [obj] = setTexture(obj, image, extrinsics, intrinsics, mask, textureDensity)
            assert(size(extrinsics,1) == 4 || size(extrinsics,2) == 4, 'Extrinsics must be a [4x4] homogeneous transformation matrix from object/model frame to camera frame' );
            assert(size(intrinsics,1) == 3 || size(intrinsics,2) == 3, 'Intrinsics must be a [3x3] camera intrinsics matrix in the standard form K = [fx 0 cx; 0 fy cy; 0 0 1]' );

            if exist('mask', 'var')
                assert(isa(mask,'logical'), 'Mask must be a two-dimensional binary array');
                obj.mask = mask;
                obj.hasMask = true;
            end
            obj.hasTexture = true;
            obj.img = image;
            obj.extrinsics = extrinsics;
            obj.intrinsics = intrinsics;
            if exist('textureDensity', 'var')
                obj.textureDensity = textureDensity;
            end
        end
        
        function [obj] = startShowingAxes(obj)
            obj.showAxis = true;
        end
        
        function [obj] = dumpOutput(obj, outputPathForPNGimage)
            obj.outputPath = outputPathForPNGimage;
            obj.hasOutputPath = true;
        end
        
        function [obj] = updatePlotPositionSize(obj, plotDetails)
            assert(numel(plotDetails) == 4, 'Plot position must be a [1 x 4] array of [px py sx sy] : (px,py) being position, (sx, sy) being [HEIGHT, WIDTH]');
            assert(size(plotDetails,1) == 1, 'Plot position must be a [1 x 4] array of [px py sx sy] : (px,py) being position, (sx, sy) being [HEIGHT, WIDTH]');
            
            obj.plotPositionAndSize = plotDetails;
        end
        
        function [T, figureHandle] = plotReconstructionNow(obj, stopNewFigure)
            
            if ~obj.hasDensePointcloud
                warning('Please set the dense pointcloud with the constructor call, returning now without doing anything');
                return;
            else
                
                fprintf('\n\n-------: Starting <strong>surface reconstruction</strong> [for vizualisation]:-------\n\n');
                P = obj.densePointCloud;

                figureHandle = [];
                
                if ~exist('stopNewFigure', 'var')
                    f1 = figure(1);
                    f1.Position = obj.plotPositionAndSize;
                    figureHandle = f1;
                    set(f1,'renderer','opengl');
                end
                
                [T]=CrustOpen(P); % actual surface reconstruction method using ball-pivot method
                
                trisurf(T,P(:,1), P(:,2), P(:,3)); hold on; shading interp; axis equal;

                if(obj.material == 0)
                    material dull;
                elseif(obj.material == 1)
                    material metal;
                elseif(obj.material == 2)
                    material shiny;                    
                end

                if(obj.hasTexture)
                    trnsfP = (obj.extrinsics(1:3,1:3)*P.' + obj.extrinsics(1:3,4)).';
                    projP = trnsfP./trnsfP(:,3);
                    [M,N] = meshgrid(linspace(min(projP(:,1)), max(projP(:,1)), obj.textureDensity), linspace(min(projP(:,2)), max(projP(:,2)), obj.textureDensity));
                    dirVec =[M(:) N(:)];
                    dirVec = [dirVec ones(size(dirVec,1),1)];
                    for ii = 1:size(dirVec,1)
                        dirVec(ii,:) = dirVec(ii,:)./norm(dirVec(ii,:));
                    end
                    tree = opcodemesh(trnsfP.',T.');
                    origin = zeros(3,size(dirVec,1));
                    [~,~,~,~,Q] = tree.intersect(origin,dirVec.');
                    imgP = ((obj.intrinsics*dirVec.').');
                    for ii = 1:size(imgP,1)
                        shouldDisplay = true;
                        px = imgP(ii,2); py = imgP(ii,1);
                        if(px<1)
                            px = 1;
                        end
                        if(py<1)
                            py = 1;
                        end
                        if(px>size(obj.mask,1))
                            px=size(obj.mask,1);
                        end
                        if(py>size(obj.mask,2))
                            py=size(obj.mask,2);
                        end
%                         fprintf('[%d, %d]\n', px, py);
                        if(obj.hasMask)
                            if~(obj.mask(round(px), round(py)))
                                shouldDisplay = false;
                            end
                        end

                        if(shouldDisplay && all(~isnan(Q(:,ii))))
                            C = interpImg(obj.img,[px,py]);
                            tQ = obj.extrinsics(1:3,1:3).'*Q(:,ii) - (obj.extrinsics(1:3,1:3).'*obj.extrinsics(1:3,4));
                            plot3(tQ(1),tQ(2),tQ(3),'o', 'MarkerFaceColor',C, 'MarkerEdgeColor','none');
                        end
                    end
                end
                
                if(obj.showAxis)
                    axis on; grid on;
                    xlabel('X', 'Interpreter', 'Latex');
                    ylabel('Y', 'Interpreter', 'Latex');
                    zlabel('Z', 'Interpreter', 'Latex');
                else
                    axis off; grid off;
                end
                
                for iL = 1:obj.numLightModel
                    light('Style',obj.lightModel{obj.numLightModel}.style,'Position',obj.lightModel{obj.numLightModel}.position); hold on;
                end
                
                if(obj.hasSparsePointcloud)
%                     pcshow(obj.sparsePointCloud, 'r', 'MarkerSize', 100); hold on;
%                     set(gcf,'color','w');
%                     set(gca,'color','w');
%                     set(gca, 'XColor', [0.0 0.0 0.0], 'YColor', [0.0 0.0 0.0], 'ZColor', [0.0 0.0 0.0]);
                    p = plot3(obj.sparsePointCloud(:,1).',obj.sparsePointCloud(:,2).',obj.sparsePointCloud(:,3).','o', 'MarkerSize', obj.sparsePointSize); hold on;
                    p.MarkerFaceColor = [0.9608 0.0941 0.4706]; p.MarkerEdgeColor = [0.9608 0.0941 0.4706];
                end
                
                if(obj.hasSparseNormals)
                    nP = size(obj.sparsePointCloud ,1);
                    
                    for iP = 1:nP
                        pC = obj.sparsePointCloud(iP, :);
                        nC = obj.sparseNormals(iP, :);
                        tC = pC + (obj.normalLength.*nC);
                        plot3([pC(1) tC(1)].', [pC(2) tC(2)].', [pC(3) tC(3)].', 'k-', 'LineWidth', 1); hold on;
                    end
                end
                
                pause(0.1); drawnow(); pause(0.1);
                
                if ( (obj.hasOutputPath) && ~exist('stopNewFigure', 'var'))
                    saveas(f1,obj.outputPath);
                end
                
                pause(0.1);
                fprintf('\n-------: Finished doing <strong>surface reconstruction</strong> [for vizualisation]:-------\n\n');
                
                return;
            end
        end
    end
    
end