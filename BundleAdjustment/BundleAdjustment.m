function [sparseReconstruction, denseReconstruction, template] = BundleAdjustment(targetPoints, Data, intrinsics, topology, nng, weights, doPlot, nDense, padding)

    nImg = size(targetPoints,2);
    nPts = size(targetPoints(1).target,1);
    
    xGridNum = 3;
    yGridNum = 3;
    zGridNum = 3;

    nSmooth = 250; nKsmooth = nSmooth-10;%round(nSmooth/6);

    if ~exist('padding', 'var')
        padding = 0.1;
    end
    shouldOutput = false;
    if exist('doPlot', 'var')
        if (doPlot == true)
            close all; 
            f = figure(1);
            f.Position = [10 20 1800 1500];
            shouldOutput = true;
        end
    end
    
    %%% Combining template

    for iF = 1:nImg
        if((topology == 0) || (topology == 1) || (topology == 2) )
            [templateN] = getNaturalCorrespondence(targetPoints(iF).template, topology);
        else
            error('Incorrect topology specified! Options are: 0 - plane, 1 - cylinder or 2 - sphere. 3 - non-parametric meshes are not handled here');
        end
        DataCell{iF} = templateN.';
        VisVecCell{iF} = ones(1,size(templateN, 1));
    end

    GPA = RigidGPA;
    GPA.addDataByArray(DataCell, VisVecCell);
    GPA.run();

    template = GPA.mPoses{1}(1:2,1:2).'*GPA.mShape - (GPA.mPoses{1}(1:2,1:2).'*GPA.mPoses{1}(:,3)) ;
    template = template.';

    %%% Re-initialising non-rigid warps
    for iF = 1:nImg
        controlGrid = getGrid(targetPoints(iF).target, xGridNum, yGridNum, zGridNum);
        nTransform = NonRigidTransformations(controlGrid);
        tpsParameterList(iF).handler = nTransform;
        controlPts(iF).target = controlGrid;
    end

    fun = @error;

    [X] = marshallParameter(controlPts, template);

    %%% Smoothing related operations
    for iF = 1:nImg
        templateInterp = getNewRandomTemplate(targetPoints(iF).template, nSmooth, topology, 0);
        smoothTemplateList{iF} = templateInterp;
    end
    [~, denseReconstruction, ~] = finaliseReconstruction(X, targetPoints, smoothTemplateList, tpsParameterList, nImg, nPts);
    for iF = 1:nImg
        DummyData.p(iF).p = denseReconstruction(iF).points.';
        [nngSmooth] = getNeighborhoodDuplicated(DummyData, nKsmooth);
        nngSmoothList{iF} = nngSmooth;
    end
    smoothData.templates = smoothTemplateList;
    smoothData.nngs = nngSmoothList;

    MAX_ITER = 5;
    TOLERANCE = 10^(-21);

    OBJECTIVE_GRADIENT = false;

%     [Jpattern] = sparsityPattern(nImg, nPts, size(nng,2));

    if(shouldOutput == true)
        options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',  ...
            'OptimalityTolerance',TOLERANCE,'StepTolerance',TOLERANCE,...
            'FunctionTolerance',TOLERANCE, 'OutputFcn',@outfun,...
            'FunValCheck','on','Diagnostics','on','MaxIterations', MAX_ITER, 'FiniteDifferenceType','forward', ...
            'CheckGradients', false, 'Display', 'iter-detailed','SpecifyObjectiveGradient', OBJECTIVE_GRADIENT);%, 'JacobPattern', Jpattern); %'SpecifyObjectiveGradient', true,
    else
        options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',  ...
            'OptimalityTolerance',TOLERANCE,'StepTolerance',TOLERANCE,...
            'FunctionTolerance',TOLERANCE, ...
            'FunValCheck','on','Diagnostics','on','MaxIterations', MAX_ITER, 'FiniteDifferenceType','forward', ...
            'CheckGradients', false, 'Display', 'iter-detailed','SpecifyObjectiveGradient', OBJECTIVE_GRADIENT);%, 'JacobPattern', Jpattern); %'SpecifyObjectiveGradient', true,
    end

    [updatedX] = lsqnonlin(fun, X, [], [], options, intrinsics, targetPoints, tpsParameterList, Data, nng, weights, topology, controlPts, smoothData);        

    if(shouldOutput == true)
        pause(1);
        close all;
    end

    for iF = 1:nImg
        templateInterp = getNewRandomTemplate(targetPoints(iF).template, nDense, topology, padding);
        templateList{iF} = templateInterp;
    end

    % Reconstruction completion
    [sparseReconstruction, denseReconstruction, template] = finaliseReconstruction(updatedX, targetPoints, templateList, tpsParameterList, nImg, nPts);
    
end

function [E] = error(X, intrinsics, targetPoints, tpsParameterList, Data, nng, weights, topology, controlPts, smoothData) 
    nImg = size(targetPoints,2);
    nPts = size(targetPoints(1).target,1);
    

    [targetUpdated, templateUpdated] = unmarshallParameters(X, nImg, nPts);

    E = []; 

    [~, denseReconstruction, ~] = finaliseReconstruction(X, targetPoints, smoothData.templates, tpsParameterList, nImg, nPts);

    for iF = 1:nImg
         nTransform = tpsParameterList(iF).handler;
         nTransform.update_tps(targetUpdated(iF).target);
         pointsInterpolated = nTransform.interpTPS(targetPoints(iF).target);

         [Ereprojection] = reprojectionError(pointsInterpolated, Data.p(iF).p, intrinsics, targetPoints(iF).T);
         [Eisometry] = isometryError(pointsInterpolated, templateUpdated, nng, topology);
         [Esmoothness] = smoothnessError(denseReconstruction(iF).points, smoothData.nngs{iF});

         Ereprojection = Ereprojection./numel(Ereprojection); 
         Eisometry = Eisometry./numel(Eisometry); 
         Esmoothness = Esmoothness./numel(Esmoothness); 

         controlPtsDiff = targetUpdated(iF).target - controlPts(iF).target;
%          controlPtsDiff = controlPtsDiff.*0;

         En = [( weights(1).*Ereprojection );( weights(2).*Eisometry );( weights(3).*Esmoothness );(1 - sum(weights))*(1/nPts)*norm(controlPtsDiff, 'fro') ];

         E = [E; En]; 
    end
end

% function [Jpattern] = sparsityPattern(nImg, nPts, kN)
%     J = zeros((nImg*(nPts + (nPts*kN) + nPts + 1)), ( (2*nPts) + (3*nImg*nPts) ));
% 
%     for iF = 1:nImg
%         rowStartIndex = (iF-1)*(nPts + (nPts*kN) + nPts) + 1;
%         rowEndIndex = rowStartIndex + nPts;
%         colStartIndex = (2*nPts) + (3*nPts*(iF-1)) + 1;
%         colEndIndex = colStartIndex + (3*nPts);
% 
%         J(rowStartIndex:rowEndIndex, colStartIndex:colEndIndex) = 1;
% 
%         rowStartIndex = rowStartIndex + nPts;
%         rowEndIndex = rowStartIndex + (nPts*kN);
% 
%         J(rowStartIndex:rowEndIndex, colStartIndex:colEndIndex) = 1;
%         J(rowStartIndex:rowEndIndex, 1:(2*nPts)) = 1;
% 
%         rowStartIndex = rowStartIndex + (nPts*kN);
%         rowEndIndex = rowStartIndex + nPts;
% 
%         J(rowStartIndex:rowEndIndex, colStartIndex:colEndIndex) = 1;
%     end
% 
%     Jpattern = J;
% end


function [stop] = outfun(X,~,~, intrinsics, targetPoints, tpsParameterList, Data, ~, ~, ~, ~, ~)
    nImg = size(targetPoints,2);
    nPts = size(targetPoints(1).target,1);

    [targetUpdated, templateUpdated] = unmarshallParameters(X, nImg, nPts);

    for iF = 1:nImg
        nTransform = tpsParameterList(iF).handler;
        nTransform.update_tps(targetUpdated(iF).target);
        points = nTransform.interpTPS(targetPoints(iF).target);

        points = targetPoints(iF).T(1:3,1:3)*points.' + targetPoints(iF).T(1:3,4);
        points = points.';

        gtPoints = Data.Pgth(iF).P;

        subplot(3,nImg,iF);
        plot3(points(:,1).',points(:,2).',points(:,3).', 'or','MarkerFaceColor','r','MarkerSize',5); hold on;
        plot3(gtPoints(1,:).',gtPoints(2,:).',gtPoints(3,:).', 'ok','MarkerFaceColor','k','MarkerSize',4);
        xlabel('X', 'Interpreter','latex','FontSize',6);
        ylabel('Y', 'Interpreter','latex','FontSize',6);
        zlabel('Z', 'Interpreter','latex','FontSize',6);
        legend({'Current reconstruction', 'GT'}, 'Interpreter', 'Latex', 'FontSize',5);
        grid on; axis on; hold off;
        pause(0.1);

        projPts = points./points(:,3);
        projPts = intrinsics*projPts.';

        subplot(3,nImg,(nImg + iF));

        scatter(projPts(1,:).',projPts(2,:).', 'r*'); hold on;
        scatter(Data.p(iF).p(1,:).',Data.p(iF).p(2,:).', 'ko'); 
        legend({'Current reprojection', 'Key-pt.s correspondence'}, 'Interpreter', 'Latex', 'FontSize',5);
        grid on; axis on; hold off;
        pause(0.1);

        subplot(3,nImg,((2*nImg) + iF));

        plot3(templateUpdated(:,1).',templateUpdated(:,2).',zeros(1, size(templateUpdated,1)), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize',4); hold on;
        plot3(points(:,1).',points(:,2).',points(:,3).', 'or', 'MarkerFaceColor', 'r', 'MarkerSize',5); hold on;
        xlabel('X', 'Interpreter','latex','FontSize',6);
        ylabel('Y', 'Interpreter','latex','FontSize',6);
        zlabel('Z', 'Interpreter','latex','FontSize',6);
        legend({'Current interp. pt.s', 'Template'}, 'Interpreter', 'Latex', 'FontSize',5);
        grid on; axis on; hold off;
        pause(0.1);

    end

    stop = 0;
    
end

function [sparseReconstruction, denseReconstruction, template] = finaliseReconstruction(updatedX, target, templateList, tpsParameterList, nImg, nPts)
    [targetUpdated, template] = unmarshallParameters(updatedX, nImg, nPts);


    for iF = 1:nImg
        nTransform = tpsParameterList(iF).handler;
        nTransform.update_tps(targetUpdated(iF).target);
        sparse = nTransform.interpTPS(target(iF).target);

        templateInterp = templateList{iF}; 
        
        nTransform = NonRigidTransformationsInterp(target(iF).template, sparse);
        dense= nTransform.interpTPS(templateInterp);        
        sparseReconstruction(iF).points = sparse;
        denseReconstruction(iF).points = dense;         
    end
end

function [template] = getNewRandomTemplate(previousAlignedPointcloud, nDense, topology, padding)
    
    [templateN] = getNaturalCorrespondence(previousAlignedPointcloud, topology);

    if(topology == 0)
        xmin = min(templateN(:,1))-padding; xmax = max(templateN(:,1))+padding;
        ymin = min(templateN(:,2))-padding; ymax = max(templateN(:,2))+padding;
        templateInterp = [random_number_within_range(xmin, xmax, nDense).' random_number_within_range(ymin, ymax, nDense).' zeros(nDense,1)];
    elseif(topology == 1)
        ymin = min(templateN(:,1))-padding; ymax = max(templateN(:,1))+padding;
        r = random_number_within_range(ymin, ymax, nDense);
        theta = random_number_within_range(-pi, pi, nDense);
        templateInterp = [sin(theta).' r.' cos(theta).'];
    elseif(topology == 2)
        theta1 = random_number_within_range(-pi, pi, nDense);
        theta2 = random_number_within_range(0, 2*pi, nDense);
        templateInterp = [(-cos(theta1).*cos(theta2)).' (-sin(theta1)).' (-cos(theta1).*sin(theta2)).'];
    end        
    template = templateInterp;
end

function [E] = reprojectionError(points, pts2Dcorrespondence, K, T)

    nPts = size(points,1);
    points = T(1:3,1:3)*points.' + T(1:3,4);
    points = points.';

    recons = points;

    points = points./points(:,3);
    points = K*points.';

    errorVec = points(1:2,:) - pts2Dcorrespondence(1:2,:);
    errorVec = errorVec.^2;

    E = sum(errorVec).';

    E = real(E);
end

function [E, J] = isometryError(points, template, nng, topology)
    nPts = size(points,1);

    EPSILON = 10^(-4);

    E = [];
    J = [];

    for iM = 1:nPts
        P1 = points(iM,:);
        T1 = template(iM,:);
        for iQ = 1:size(nng,2)
            P2 = points(nng(iM,iQ),:);
            T2 = template(nng(iM,iQ),:);

            if(topology == 2)
                [geodesic] = arcLength(T1, T2); % for sphere, Delta is not isometric
            else
                geodesic = EuclideanDistance(T1,T2);
            end

            tau_1 = sqrt( (P1(1) - P2(1))^2 + (P1(2) - P2(2))^2 + (P1(3) - P2(3))^2  + EPSILON);
            tau_2 = sqrt(geodesic^2 + EPSILON);

            error = real(tau_1 - tau_2);


            E = [E; error];
        end
    end

    E = real(E);
end

function [E] = smoothnessError(points, nng)
    nPts = size(points,1);
    % Following Gaussian smoothing: https://graphics.stanford.edu/courses/cs468-01-fall/Papers/taubin-smoothing.pdf
    % weights are our custom implementation
    E = [];
    for iM = 1:nPts
        P1 = points(iM,:);
        error = 0;
        for iQ = 1:(size(nng,2))
            P2 = points(nng(iM,iQ),:);
            V1 = full(P2 - P1);
            error = error + (V1); % 1/d makes the cost 1
        end
        E = [E; (error).'];
    end
    E = real(E);
end

function [X] = marshallParameter(target, template)
    nImg = size(target,2);
    nPts = size(target(1).target,1);
    Y1 = [];
    for iF = 1:nImg
        Y1 = [Y1;target(iF).target];
    end
    X = [reshape(template,[2*size(template,1),1]);reshape(Y1,[3*nPts*nImg,1])];
end

function [target, template] = unmarshallParameters(X, nImg, nPts)
    lenTempl = nPts*2;
    lenGrid = (numel(X) - lenTempl)/(3*nImg);
    Y = reshape(X(lenTempl+1:end),[(nImg*lenGrid), 3]);
    template = reshape(X(1:lenTempl),[nPts,2]);
    for iF = 1:nImg
        startIndex = (lenGrid*(iF-1)) + 1;
        endIndex = lenGrid*iF;
        target(iF).target = Y(startIndex:endIndex,:);
    end
    
end

