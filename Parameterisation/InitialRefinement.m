
function [interpolated_data, T, cleaned_data, warpInfo]...
    = InitialRefinement(parameterised_reconstruction_target, visibility, max_iter, topology, NA_Model,... 
                                            alpha_param, beta_param, interpolation_density, padding, debugDisplay, cylinderReferenceVector)

    fprintf('\n\n-------: Starting <strong>initial parameterization</strong>:-------\n');
    global tps_g N_target_g Linv_g xi_g base_g aligned_reconstruction_g alpha beta topology_g non_analytic_target;

    if ~exist('debugDisplay', 'var')
        debugDisplay = false;
    end

    if(topology == 0)
        tps_g = WPT;
        align = PlaneAlignment;
    elseif(topology == 1)
        tps_g = WTT;
        align = CylinderAlignment;
    elseif(topology == 2)
        tps_g = WST;
        align = SphereAlignment;
    else
        error('Wrong topology specified!');
    end

    if ~exist('padding', 'var')
        padding = 0.1;
    end    

    if ~exist('interpolation_density', 'var')
        interpolation_density = 5000;
    end
    topology_g = topology;

    alpha = alpha_param;
    beta = beta_param;

    [parameterised_reconstruction_target, discardedList] = remove_invalid_points(parameterised_reconstruction_target, visibility);

    if( (topology == 1) && exist('cylinderReferenceVector', 'var'))
        [aligned_reconstruction_g, r, C, T] = align.align_pointcloud(parameterised_reconstruction_target, cylinderReferenceVector);
    else
        [aligned_reconstruction_g, r, C, T] = align.align_pointcloud(parameterised_reconstruction_target);
    end
    if(topology == 3)
        non_analytic_target = aligned_reconstruction_g;
        base_g = align.flattened_target;
        aligned_points = aligned_reconstruction_g(:,1:3);
        data = aligned_points(:);
    else
        [base_g] = tps_g.find_correspondence_natural(aligned_reconstruction_g, length_matrix(aligned_reconstruction_g), 1);
%         data = aligned_reconstruction_g(:); % for cylinders (at least) -> use this!
        data = base_g(:);
    end



    template = base_g;
    N_target_g = length_matrix(aligned_reconstruction_g);


    [A, d, Linv, xi] = tps_g.compute_tps(base_g, aligned_reconstruction_g, N_target_g);
    Linv_g = Linv; xi_g = xi;

    fun = @error;

    if(debugDisplay)
        close all;
        figure;
        options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',  ...
            'FiniteDifferenceStepSize',0.001,'StepTolerance',10.^(-12),  ...
            'FunctionTolerance',10.^(-12),'Display', 'iter-detailed', 'OutputFcn',@outfun,  ...
            'MaxIterations',max_iter,'FunValCheck','off','Diagnostics','off',  ...
            'SpecifyObjectiveGradient', true,'UseParallel',false);
    else
        options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',  ...
            'FiniteDifferenceStepSize',0.001,'StepTolerance',10.^(-12),  ...
            'FunctionTolerance',10.^(-12),'Display', 'iter-detailed',  ...
            'MaxIterations',max_iter,'FunValCheck','off','Diagnostics','off',  ...
            'SpecifyObjectiveGradient', true,'UseParallel',false);
    end

    lb = []; ub = [];

    [X, resnorm, residual, exitflag, output] = lsqnonlin(fun, data, lb, ub, options);
    [interpolated_data_sparse, E, J, jacobian_components] = interpolate_reconstruction(X);
    cleaned_data = interpolated_data_sparse;
    fprintf('\n\n======= Starting <strong>high-resolution interpolation</strong>, this may take a few minutes =======\n\n');
    [interpolated_data, E, denseTemplate] = densely_interpolate_reconstruction(X, interpolation_density, topology, padding, NA_Model);
    fprintf('\n\n======= Done <strong>high-resolution interpolation</strong> =======\n\n');

    if(debugDisplay)
        close all;
    end

    fprintf('\n-------: Finished doing <strong>initial parameterization</strong>:-------\n\n');

    warpInfo.discardedList = discardedList;
    warpInfo.xi = xi;
    warpInfo.Linv = Linv;
    warpInfo.template = template;
end

function [E,J] = error(data)

    global N_target_g base_g aligned_reconstruction_g alpha beta non_analytic_target topology_g;

    [interpolated_data, E, J, jacobian_components] = interpolate_reconstruction(data);

    if(topology_g == 3)
        [E_eucl, J_eucl] = error_point_to_plane(interpolated_data, non_analytic_target, jacobian_components, N_target_g);
        [E_norm, J_norm] = error_norm_bending(interpolated_data, aligned_reconstruction_g, jacobian_components, N_target_g);                
        [E_eucl_b, J_eucl_b] = error_point_to_plane(interpolated_data, aligned_reconstruction_g, jacobian_components, N_target_g);
    else
        [E_eucl, J_eucl] = error_point_to_plane(interpolated_data, aligned_reconstruction_g, jacobian_components, N_target_g);
        [E_norm, J_norm] = error_norm_bending(interpolated_data, base_g, jacobian_components, N_target_g);        
    end

    
    N_eucl = size(E_eucl,1); N_norm = size(E_norm,1); N_bend = size(E,1);

    if(topology_g == 3)
        N_bend = size(E_eucl_b,1);
        E = vertcat(alpha*E_eucl*(1/N_eucl), beta*E_norm*(1/N_norm), (0)*E_eucl_b*(1/N_bend)); 
        J = vertcat(alpha*J_eucl*(1/N_eucl), beta*J_norm*(1/N_norm), (0)*J_eucl_b*(1/N_bend));
    else
        E = vertcat(alpha*E_eucl*(1/N_eucl), beta*E_norm*(1/N_norm), (1 - alpha - beta)*E*(1/N_bend)); 
        J = vertcat(alpha*J_eucl*(1/N_eucl), beta*J_norm*(1/N_norm), (1 - alpha - beta)*J*(1/N_bend));
    end
end

function [stop] = outfun(data, ~, ~)
    global aligned_reconstruction_g;
    [interpolated_data, ~, ~, ~] = interpolate_reconstruction(data);
    points = full(interpolated_data(:,1:3));
    normals = full(interpolated_data(:,4:6));

    plot3(points(:,1).',points(:,2).',points(:,3).','ro'); hold on;
    plot3(aligned_reconstruction_g(:,1).',aligned_reconstruction_g(:,2).',aligned_reconstruction_g(:,3).','bo'); hold on;
    for iP = 1:size(normals,1)
        P1 = points(iP,:);
        P2 = P1 + (0.2.*normals(iP,:));
        plot3([P1(1) P2(1)].',[P1(2) P2(2)].',[P1(3) P2(3)].','k-', 'HandleVisibility','off'); hold on;
    end
    legend('Interp. pt.s', 'Aligned recons.');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on; axis on;
    hold off;
    pause(0.01); drawnow();
    stop = 0;
end

function [E,J] =  error_point_to_plane(interpolated_data, aligned_reconstruction, jacobian_components, N_target)    
    E = zeros(N_target,1);
    J = zeros(N_target,(3*N_target));

    for index = 1: N_target
        E(index,1) =      ((interpolated_data(index,1)-aligned_reconstruction(index,1))*interpolated_data(index,4) + ...
            (interpolated_data(index,2)-aligned_reconstruction(index,2))*interpolated_data(index,5) + ...
            (interpolated_data(index,3)-aligned_reconstruction(index,3))*interpolated_data(index,6) );

        mu_j = [(interpolated_data(index,1)-aligned_reconstruction(index,1)) (interpolated_data(index,2)-aligned_reconstruction(index,2)) (interpolated_data(index,3)-aligned_reconstruction(index,3)) 0];
        J1 = jacobian_components.first_term(index).data;
        dN_dY = jacobian_components.second_term(index).data;
        J_row = J1 + (mu_j*dN_dY);
        J(index,:) = J_row;
    end

end

function [E,J] =   error_norm_bending(interpolated_data, base, jacobian_components, N_target)
    global topology_g concavity_g;
    E = zeros(N_target,1);
    J = zeros(N_target,(3*N_target));

    for index = 1: N_target
        if(topology_g == 0)
            N_base = ([0 0 1] - base(index,:))/norm(base(index,:));
        else
            if(topology_g == 1)
                base(index,2) = 0.0;
            end
            N_base = -base(index,:)/norm(base(index,:));
        end
        nx_ = -N_base(1,1) + interpolated_data(index,4) ;
        ny_ = -N_base(1,2) + interpolated_data(index,5) ;
        nz_ = -N_base(1,3) + interpolated_data(index,6) ;
        
        E(index,1) =  sqrt(nx_.^2 + ny_.^2 + nz_.^2);

        dN_dY = jacobian_components.second_term(index).data;

        dN_dY_x = dN_dY(1,:);
        dN_dY_y = dN_dY(2,:);
        dN_dY_z = dN_dY(3,:);

        J_row =  (dN_dY_x * nx_) + (dN_dY_y * ny_) + (dN_dY_z * nz_);
        J(index,:) = J_row./E(index,1);
    end
end


function [interpolated_data, E, J, jacobian_components] = interpolate_reconstruction(control_handles)
    global tps_g N_target_g Linv_g xi_g base_g;
    target = reshape(control_handles, [N_target_g 3]);
    [A, d, E,J] = tps_g.update_tps(target, Linv_g, xi_g, N_target_g);
    [X, Y, Z, Nx, Ny, Nz, E_temp, jacobian_components] = tps_g.interpolate_tps(A, d, ...
        base_g, base_g, N_target_g, N_target_g, xi_g, Linv_g, 0);
    interpolated_data = horzcat(X, Y, Z, Nx, Ny, Nz);
end

function [interpolated_data, E, targetInterp] = densely_interpolate_reconstruction(control_handles, numPts, topology, padding, NA_Model)
    global tps_g N_target_g Linv_g xi_g base_g;
    target = reshape(control_handles, [N_target_g 3]);
    if (topology == 0)
        xmin = min(target(:,1))-padding; xmax = max(target(:,1))+padding;
        ymin = min(target(:,2))-padding; ymax = max(target(:,2))+padding;
        targetInterp = [random_number_within_range(xmin, xmax, numPts).' random_number_within_range(ymin, ymax, numPts).' zeros(numPts,1)];
    elseif (topology == 1)
        ymin = min(target(:,2))-padding; ymax = max(target(:,2))+padding;
        r = random_number_within_range(ymin, ymax, numPts);
        theta = random_number_within_range(-pi, pi, numPts);
        targetInterp = [sin(theta).' r.' cos(theta).'];
    elseif (topology == 2)
        theta1 = random_number_within_range(-pi, pi, numPts);
        theta2 = random_number_within_range(0, 2*pi, numPts);
        targetInterp = [(-cos(theta1).*cos(theta2)).' (-sin(theta1)).' (-cos(theta1).*sin(theta2)).'];
    else
        xmin = min(base_g(:,1)); xmax = max(base_g(:,1));
        ymin = min(base_g(:,2)); ymax = max(base_g(:,2));
        targetInterp = [random_number_within_range(xmin, xmax, numPts).' random_number_within_range(ymin, ymax, numPts).' zeros(numPts,1)];
    end
    
    [A, d, E,J] = tps_g.update_tps(target, Linv_g, xi_g, N_target_g);
    [X, Y, Z, Nx, Ny, Nz, ~, ~] = tps_g.interpolate_tps(A, d, ...
        targetInterp, base_g, numPts, N_target_g, xi_g, Linv_g, -2);
    interpolated_data = horzcat(X, Y, Z, Nx, Ny, Nz);
end

function [N] = length_matrix(X)
    N = size(X);
    N = N(1);
end

function [data_cleaned, discarded_indices_retVal] = remove_invalid_points(data, visibility)
%     [data_cleaned_outlier, inlier_indices, discarded_indices_retVal]  = pcdenoise(pointCloud(data));
    import java.util.*;
    N = size(data);
    N = N(1);
    discarded_indices = Stack();

    for index = 1:N
        if(isnan(data(index,1)) || isnan(data(index,2)) || isnan(data(index,3)) || (visibility(1,index) == 0)  ) %|| ismember(index,discarded_indices_retVal) 
            discarded_indices.push(index);
        end
    end

    S = size(discarded_indices);
    discarded_indices_retVal = zeros(S,1);

    for index = 1:S
        p = discarded_indices.pop();
        data(p,:) = [];
        discarded_indices_retVal(index,1) = p;
    end
    data_cleaned = data;   
end

function [NA_Model] = adjustModelScale(reconstructedPointcloud, NA_Model)
    diagonalModel = getLargestDiagonal(NA_Model.model.Location);
    diagonalFlatModel = getLargestDiagonal(NA_Model.flat_model.Location);
    diagonalPointcloud = getLargestDiagonal(reconstructedPointcloud);

    scaleModel = diagonalPointcloud/diagonalModel;
    scaleFlatModel = diagonalPointcloud/diagonalFlatModel;

    scaledModel = scaleModel.*NA_Model.model.Location;

    minZ = min(scaledModel(:,3));

    if(minZ < 0)
        scaledModel = scaledModel - [0 0 (minZ - 0.01)];
    end

    NA_Model.model = pointCloud(scaledModel);
    NA_Model.flat_model = pointCloud(scaleFlatModel.*NA_Model.flat_model.Location);
end
