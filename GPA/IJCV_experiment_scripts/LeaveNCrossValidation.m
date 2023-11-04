function predictedPointClouds = LeaveNCrossValidation (DataCell, VisVecCell, dimFeatureSpace, smoothParam, N, gaugeShape)

nShapes =  size(DataCell, 2);

mPoints = size(DataCell{1}, 2);

dDim = size(DataCell{1}, 1);

predictedPointClouds = DataCell; % initialize storage. same size as input datum shapes

nTrials = ceil(mPoints / N);

for kk = 1 : nTrials
    
    AllPointsFlag = ones(1, mPoints);
    
    testingFragment = (1 + (kk-1) * N)  : min(kk*N, mPoints);
    
    AllPointsFlag (testingFragment) = 0;
    
    trainingPointsFlag = logical(AllPointsFlag);
    
    testingPointsFlag = logical(ones(1, mPoints) - AllPointsFlag);
    
    TestPointsPrediction = runTrainTestData (DataCell, VisVecCell, dimFeatureSpace, smoothParam, trainingPointsFlag, testingPointsFlag, gaugeShape);
    
    for ii = 1 : nShapes
        
        predictedPointClouds{ii}(:, testingPointsFlag) =  TestPointsPrediction{ii};
        
    end
    
end

end





function TestPointsPrediction = runTrainTestData (DataCell, VisVecCell, dimFeatureSpace, smoothParam, trainingPointsFlag, testingPointsFlag, gaugeShape)

    if (smoothParam == inf && dimFeatureSpace == 0)
        GPA = DefGPA('AFFINE');
    elseif (smoothParam > 0 && dimFeatureSpace > 0)
        GPA = DefGPA('TPS');
    elseif (smoothParam == 0  && dimFeatureSpace == -1)
        GPA = GPA_Functor('EUC-ALL');
    elseif (smoothParam == -1  && dimFeatureSpace == -1)
        GPA = GPA_Functor('AFF-ALL');
    else
        fprintf(2, '\nwrong smoothParma and controlPointsPerAxis, cannot identify a proper algorithm to test\n');
    end

    nShapes =  size(DataCell, 2);
    
    TestPointsPrediction = cell(1, nShapes);
    
    for ii = 1 : nShapes
        
        GPA.addPointCloud (DataCell{ii}(:, trainingPointsFlag), VisVecCell{ii}(:, trainingPointsFlag));
        
    end
    
    GPA.dimFeatureSpace = dimFeatureSpace;    

    GPA.flagJointlyOptimizeLamdaRigidTransformation = false;

    GPA.smoothParam = smoothParam;
     
    GPA.run(smoothParam);
  
    [R, t] = EuclideanProcrustes (GPA.mShape, gaugeShape(:, trainingPointsFlag));
    
    for ii = 1 : nShapes
        
        TestPointsPrediction{ii} = R * GPA.transformPoints(DataCell{ii}(:, testingPointsFlag), ii) + t;
        
        TestPointsPrediction{ii} = TestPointsPrediction{ii} .* VisVecCell{ii}(:, testingPointsFlag);
        
    end

    clear GPA;

end




% The solution to the similarity Procrustes problem
% sR, t = minimize || sR * D1  + t  - D2 ||
function [sR, t] = SimilairtyProcrustes (D1, D2)

meanD1 = mean(D1, 2);

meanD2 = mean(D2, 2);

M = (D1 - meanD1) * (D2 - meanD2)';

N = (D1 - meanD1) * (D1 - meanD1)';

[U, ~, V] = svd(M);

R = V * U';

if det(R) < 0
    
    fprintf(2, '\nThere exists reflection between solutions! \n');
    
end

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

if det(R) < 0
    
    fprintf(2, '\nThere exists reflection between solutions! \n');
    
end

t =  meanD2 - R * meanD1;

end
