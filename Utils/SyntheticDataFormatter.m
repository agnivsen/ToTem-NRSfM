classdef SyntheticDataFormatter<handle

    properties(Access = private)
        doShowPlots = false;
        maxImages = 6;
        maxTrainFeats = 50;
        maxTestFeats = 100;
        addedNoiseInPixels = 0.0;
        dataPath = '';
        dataHasBeenSet = false;
        T = eye(4);
    end

    methods(Access = public)
        function [obj] = SyntheticDataFormatter(dataFolderPath)
            obj.dataPath = dataFolderPath;
            obj.dataHasBeenSet = true;
        end

        function [obj] = setMaxImages(obj, maxImg)
            obj.maxImages = maxImg;
        end

        function [obj] = setMaxTrainingFeatures(obj, maxFeat)
            obj.maxTrainFeats = maxFeat;
        end

        function [obj] = setMaxTestFeatures(obj, maxFeat)
            obj.maxTestFeats = maxFeat;
        end

        function [obj] = addNoiseToFeaturesInPixels(obj, noise)
            obj.addedNoiseInPixels = noise;
        end

        function [obj] = startPlotting(obj)
            obj.doShowPlots = true;
        end

        function [obj] = stopPlotting(obj)
            obj.doShowPlots = false;
        end

        function [obj] = setGlobalTransform(obj, globalT)
            obj.T = globalT;
        end

        function [TrainData, TestData, intrinsics, transforms] = simulateSyntheticData(obj)
            if ~obj.dataHasBeenSet
                error('Call the constructor of the <SyntheticDataFormatter> class and set the data folder path! Cannot proceed otherwise.');
            end

            numFiles = dir([obj.dataPath '/*.stl']);
            numFiles = size(numFiles,1);

            assert((numFiles > 0) && (mod(numFiles,2) == 0),'There needs to be *.stl files in the input data folder, each image/model paired as <i_test.stl> and <i_train.stl> for the i-th frame');

            numFiles = numFiles/2;

            if(numFiles > obj.maxImages)
                imgIndices = randperm(numFiles);
                imgIndices = imgIndices(1:obj.maxImages);
                numFiles = obj.maxImages;
            else
                imgIndices = 1:numFiles;
            end

            f = random_number_within_range(600,800,1);

            intrinsics = [f 0 960; 0 f 540;0 0 1];

            trainPath = strcat(obj.dataPath, num2str(imgIndices(1)), '_train.stl');
            testPath = strcat(obj.dataPath, num2str(imgIndices(1)), '_test.stl');

            [trainingVertices, trainingFaces, ~] = stlReadFileEx(trainPath);
            [testingVertices, testingFaces,~] = stlReadFileEx(testPath);

            numTrainingFeats = size(trainingVertices,1);
            numTestingFeats = size(testingVertices,1);

            numTrainingOrg = numTrainingFeats;
            numTestOrg = numTestingFeats;

            if(numTrainingFeats > obj.maxTrainFeats)
                indicesTrain = randperm(numTrainingFeats);
                indicesTrain = indicesTrain(1:obj.maxTrainFeats);
                numTrainingFeats = obj.maxTrainFeats;
            else
                indicesTrain = 1:numTrainingFeats;
            end

            if(numTestingFeats > obj.maxTestFeats)
                indicesTest = randperm(numTestingFeats);
                indicesTest = indicesTest(1:obj.maxTestFeats);
                numTestingFeats = obj.maxTestFeats;
            else
                indicesTest = 1:numTestingFeats;
            end

            fprintf('\n\n\n===============================================================================\n');
            fprintf('Simulating data from Blender generated <strong>synthetic models</strong> with the following parameters:\n');
            fprintf('\t1)\tNo. of <strong>images</strong> = %d\n', numFiles);
            fprintf('\t2)\tNo. of training <strong>features</strong> = %d\n', numTrainingFeats);
            fprintf('\t3)\tNo. of testing <strong>features</strong> = %d\n', numTestingFeats);
            fprintf('\t4)\t<strong>Focal length</strong> of virtual camera = %d\n', f);
            if(obj.addedNoiseInPixels > 0.00000001)
                fprintf('\t5)\tAdded <strong>noise</strong> (in px.) = %d\n', obj.addedNoiseInPixels );
            end
            fprintf('===============================================================================\n\n\n');

            if(obj.doShowPlots)
                f1 = figure;
            end

            for iF = 1:numFiles
                trainPath = strcat(obj.dataPath, num2str(imgIndices(iF)), '_train.stl');
                testPath = strcat(obj.dataPath, num2str(imgIndices(iF)), '_test.stl');

                [trainingVertices, trainingFaces, ~] = stlReadFileEx(trainPath);
                [testingVertices, testingFaces,~] = stlReadFileEx(testPath);

                numTrainingFeatsI = size(trainingVertices,1);
                numTestingFeatsI = size(testingVertices,1);

                assert(numTrainingFeatsI == numTrainingOrg,'No. of vertices of the input STL files cannot change!');
                assert(numTestingFeatsI == numTestOrg,'No. of vertices of the input STL files cannot change!');

                trainingVertices = trainingVertices(indicesTrain,:);
                testingVertices = testingVertices(indicesTest,:);

                R = rotx(random_number_within_range(0,0.26,1))*roty(random_number_within_range(0,0.26,1))*rotz(random_number_within_range(1,1.04,1));
                t = random_number_within_range(0,1,3).' + [0;0;random_number_within_range(3,5,1)];

                T2 = [R t; 0 0 0 1];

                T = T2*obj.T;

                R = T(1:3,1:3); t = T(1:3,4);

                transforms{iF} = T;

                trainingVertices = R*trainingVertices.' + t;
                trainingVertices = trainingVertices.';

                testingVertices = R*testingVertices.' + t;
                testingVertices = testingVertices.';


                trainingVerticesProj = trainingVertices./trainingVertices(:,3);
                trainingVerticesProj = intrinsics*trainingVerticesProj.';

                trainingVerticesProj = trainingVerticesProj(1:2,:) + (obj.addedNoiseInPixels .* rand(2,numTrainingFeats));
                trainingVerticesProj = [trainingVerticesProj; ones(1, numTrainingFeats)];

                TrainData.Pgth(iF).P = trainingVertices.';
                TrainData.p(iF).p = trainingVerticesProj;

                if(obj.doShowPlots)
                    subplot(2,numFiles,iF);
                    plot3(testingVertices(:,1).',testingVertices(:,2).',testingVertices(:,3).', 'ok','MarkerFaceColor','k','MarkerSize',4); hold on;
                    plot3(trainingVertices(:,1).',trainingVertices(:,2).',trainingVertices(:,3).', 'or','MarkerFaceColor','r','MarkerSize',6); hold on;
                    xlabel('X', 'Interpreter','latex','FontSize',6);
                    ylabel('Y', 'Interpreter','latex','FontSize',6);
                    zlabel('Z', 'Interpreter','latex','FontSize',6);
                    grid on; axis on;
                    hold off; pause(0.1);

                    subplot(2,numFiles,(numFiles + iF));
                    scatter(trainingVerticesProj(1,:), trainingVerticesProj(2,:), 'r*');
                    pause(0.1);
                end

                TestData.Pgth(iF).P = testingVertices.';
            end

            TrainData.v = ones(numFiles, numTrainingFeats);

            if(obj.doShowPlots)
                pause(0.5);
                close(f1);
            end

        end


    end
end