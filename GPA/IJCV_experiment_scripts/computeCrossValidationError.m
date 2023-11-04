function cve = computeCrossValidationError(DataCell, VisVecCell, dimFeatureSpace, smoothParam, referenceShape)

cve = 0;

if(size(DataCell{1}, 2) < 500)
    testsetSize = 1;
else
    testsetSize = 40;
end

predictedPointClouds = LeaveNCrossValidation (DataCell, VisVecCell, dimFeatureSpace, smoothParam, testsetSize, referenceShape);

n = size(predictedPointClouds, 2);
m = size(referenceShape, 2);

normalizer = 0;

for ii = 1 : n
    
    pointsDeviation = (predictedPointClouds{ii} - referenceShape) .* VisVecCell{ii};
    
    sqrtError =  norm(pointsDeviation, 'fro');
    
    cve = cve + sqrtError * sqrtError;
    
    normalizer = normalizer + nnz(VisVecCell{ii});
    
end

cve = sqrt(cve/(normalizer));

end
