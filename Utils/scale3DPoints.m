function [scaledPoints, delIndex] = scale3DPoints(pointsToScale, scaleReferencePoints, doTranslate)
    %% scale3DPoints scale 3D points to a reference scale, ignoring NaN, 0, Inf points
     % *PARAMETERS*
     % * _pointsToScale_: [N x 3] in R^3  (pointcloud A)
     % * _scaleReferencePoints_: [N x 3] in R^3  (pointcloud B)
     % * _doTranslate_: boolean variable indicating if the pointcloud needs to be translated while registering with ABSOR
     %
     % Pointcloud A is scaled to match the scale of pointcloud B

     assert(size(pointsToScale,2) == 3,'scale3DPoints: input arguments should be [N x 3] points');
     assert(size(scaleReferencePoints,2) == 3,'scale3DPoints: input arguments should be [N x 3] points');
     assert(size(scaleReferencePoints,1) == size(pointsToScale,1),'scale3DPoints: input arguments should be [N x 3] points and of same dimension');
     
     shouldTranslate = false;
     if exist('doTranslate','var')
         if(doTranslate)
            shouldTranslate = true;
         end
    end
     
    Pest = pointsToScale;
    Pgt = scaleReferencePoints;

    delIndex = [];
    toDelete = false;

    for iP = 1:size(Pgt,1)
        if( (abs(Pgt(iP, 3)) < 10^(-20)) || isnan(Pgt(iP, 3)) || isinf(Pgt(iP, 3)) ||...
                isnan(Pgt(iP, 2)) || isinf(Pgt(iP, 2)) || isnan(Pgt(iP, 1)) || isinf(Pgt(iP, 1)) )
            delIndex = [delIndex iP];
            toDelete = true;
        end
    end

    if(toDelete)
        Pgt(delIndex,:) = [];
        Pest(delIndex,:) = [];
    end
    
    if(shouldTranslate)
        [regParams] = absor(Pest.', Pgt.','doScale',true,'doTrans', true);
        scaledPoints = regParams.s*regParams.R*pointsToScale.' + regParams.t;    
        scaledPoints = scaledPoints.';
    else
        [regParams] = absor(Pest.', Pgt.','doScale',true,'doTrans', false);
        scaledPoints = regParams.s.*pointsToScale;
    end
end