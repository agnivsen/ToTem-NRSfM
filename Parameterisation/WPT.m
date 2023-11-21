classdef WPT <handle
    methods
        
        %%
        %%%%%%%%%%%%%%%%%%%%
        % Projects a set of 3D points (aligned to the origin of R3) to the
        % surface of a plane
        %
        % @param: points - [Nx3] 3D points
        % @param: N - size of points
        % @param: radius - unused variable for planar (has been used in WTT)
        % @return: [corresponded_points] - points lying on surface of plane
        %%%%%%%%%%%%%%%%%%%%
        function [corresponded_points] = find_correspondence_natural(obj, points, N, radius)
            
            corresponded_points = zeros(N,3,'like',speye(2));
            
            for index = 1:N
                corresponded_points(index,1) = points(index,1);
                corresponded_points(index,2) = points(index,2);
                corresponded_points(index,3) = 0.0;
            end
        end
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%
        %  Find the (x, y) coordinates of a point lying on a plane
        %
        % @param: X, Y, Z - cartesian coordinates
        %
        % @return: [phi, theta] - computed values (here phi = x, theta = y, the names remain unchanged from WTP, to allow easier porting of code)
        %%%%%%%%%%%%%%%%%%%%
        function [phi, theta] = flattened_coordinates(obj, X, Y, Z)
            if(abs(Z) < 0.0001)
                phi = X;
                theta = Y;
            else
                warning('Asked to flatten 3D point which do not exactly lie on the XY-plane.');
                fprintf('\t---> This method is best avoided for points not in XY-plane!\n');
            end
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%
        %  Find the cartesian coordinates of a point from its 2D flattened
        %  coordinates
        %
        % @param: x, y - coordinates on flat template
        %
        % @return: [X, Y, Z] - computed cartesian coordinates
        %%%%%%%%%%%%%%%%%%%%
        function [X, Y, Z] = cartesian_coordinates(obj, x, y)
            X = x;
            Y = y;
            Z = 0;
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%
        %  First derivative of the canonical warp (Delta) w.r.t coordinates
        %  of the 2D template
        %
        % @param: x, y - spherical parameterisation
        %
        % @return: [delta_theta_1, delta_theta_2] - computed derivatives, i.e.,
        %                                                                   2 sets of [1x3] vectors
        %%%%%%%%%%%%%%%%%%%%
        function [delta_x, delta_y] = cartesian_coordinates_first_derivative(obj, x, y)
            delta_x = [1 0 0 1];  %Changes needed in Line: 603 as well
            delta_y = [0 1 0 1];
        end    

        %%
        %%%%%%%%%%%%%%%%%%%%
        %  Interpolate the 3D points from undeformed, unit sphere
        %  based on previously obtained TPS values
        %
        % @param: A, d - TPS parameter
        % @param: interpolation_candidate_vertices - vertices on unit sphere
        %         that needs to be interpolated
        % @param: control_handles - undeformed/original control handles defined
        %                           for the TPS
        % @param: N_test - number of vertices in the test data (which will be
        %                     interpolated)
        % @param: N_control - number of 3D points in data
        %
        % @param: Linv - from TPS computation
        %
        % @param: mode - used to control the type of Jacobians computed...
        %                             0: initial refinement related
        %                             Jacobian components
        %                             1: global refinement related Jacobian
        %                             components
        %
        % @return: [Px, Py, Pz, Nx, Ny, Nz] - interpolated points & normals
        %%%%%%%%%%%%%%%%%%%%
        function [Px, Py, Pz, Nx, Ny, Nz, E, jacobian_components, dQ_dxi, C, dC_dXi] = interpolate_tps(obj, A, d, interpolation_candidate_vertices, ...
                control_handles, N_test, N_control, xi, Linv, mode)


            % mode: -1 : just points, nothing else
            % mode: -2: just points & normals, nothing else
            if(mode<=-1)
                Px =  zeros(N_test, 1);
                Py =  zeros(N_test, 1);
                Pz =  zeros(N_test, 1);
                Nx = zeros(N_test, 1);
                Ny = zeros(N_test, 1);
                Nz = zeros(N_test, 1);
                for i = 1:N_test
                    P_t = ones(1, 4, 'double');
                    P_t(1,1) = interpolation_candidate_vertices(i,1);
                    P_t(1,2) = interpolation_candidate_vertices(i,2);
                    P_t(1,3) = interpolation_candidate_vertices(i,3);
                    Phi_t = zeros(1, N_control);

                    for j = 1:N_control
                        Phi_t(1,j) = obj.U(  control_handles(j,1) -  P_t(1,1),  ...
                            control_handles(j,2) -  P_t(1,2), ...
                            control_handles(j,3) -  P_t(1,3)  );
                    end
                    fp = (P_t*d) + (Phi_t*A);
                    if(mode == -2)
                        [r, theta] = obj.flattened_coordinates(P_t(1,1), P_t(1,2), P_t(1,3));
    
                        [Nx(i,1), Ny(i,1), Nz(i,1), ~, ~, ~, ~, ~, ~, ~] = obj.get_normal(A, d, control_handles, r, theta, P_t(1,1), P_t(1,2), P_t(1,3), N_control);
                    end

                    Px(i,1) = fp(1,1);
                    Py(i,1) = fp(1,2);
                    Pz(i,1) = fp(1,3);
                end

                E = [];
                jacobian_components = [];
                dQ_dxi = [];
                C = [];
                dC_dXi = [];

                return;
            end

            Px =  zeros(N_test, 1, 'like',speye(2));
            Py =  zeros(N_test, 1, 'like',speye(2));
            Pz =  zeros(N_test, 1, 'like',speye(2));
            Nx = zeros(N_test, 1, 'like',speye(2));
            Ny = zeros(N_test, 1, 'like',speye(2));
            Nz = zeros(N_test, 1, 'like',speye(2));
            E = zeros(N_test*3, 1, 'like',speye(2));
            C = zeros(N_test,1, 'like',speye(2));
            dQ_dxi = zeros(N_test, (12*N_test),'like',speye(2));
            dC_dXi = zeros(N_test, (3*N_control),'like',speye(2));
            
            dQ_dxi = sparse(dQ_dxi);
            dC_dXi = sparse(dC_dXi);
            
            if(mode == 1) % primarily used for global refinement
                

                dYT_dxi = zeros(4,(3*N_control),'like',speye(2));
                offset = 0;
                
                for loop = 1:(N_control)
                    dYT_dxi(1,loop+(3*offset*N_control)) = 1.0;
                    dYT_dxi(2,loop+(3*offset*N_control) + N_control) = 1.0;
                    dYT_dxi(3,loop+(3*offset*N_control) + N_control + N_control) = 1.0;
                    offset = offset + 1;
                end
                
                I_3l = speye(3*N_control);
                d_phi_rrT_dxi_part1 = dYT_dxi*(kron(xi, I_3l));
            end
            
            delta_theta1_theta1 = zeros(1,3,'like',speye(2));
            delta_theta2_theta2 = zeros(1,3,'like',speye(2));
            delta_theta1_theta2 = zeros(1,3,'like',speye(2));
            
            for i = 1:N_test
                P_t = ones(1, 4, 'double');
                P_t(1,1) = interpolation_candidate_vertices(i,1);
                P_t(1,2) = interpolation_candidate_vertices(i,2);
                P_t(1,3) = interpolation_candidate_vertices(i,3);
                Phi_t = zeros(1, N_control, 'like',speye(2));
                
                for j = 1:N_control
                    Phi_t(1,j) = obj.U(  control_handles(j,1) -  P_t(1,1),  ...
                        control_handles(j,2) -  P_t(1,2), ...
                        control_handles(j,3) -  P_t(1,3)  );
                end
                fp = (P_t*d) + (Phi_t*A);
                
                [r, theta] = obj.flattened_coordinates(P_t(1,1), P_t(1,2), P_t(1,3));
                
                [Nx(i,1), Ny(i,1), Nz(i,1), lambda_r_prio, lambda_theta_prio, phi_r, phi_theta, N_norm, cT, sT] = obj.get_normal(A, d, control_handles, r, theta, P_t(1,1), P_t(1,2), P_t(1,3), N_control);
                lambda_r_prio = sparse(lambda_r_prio);
                lambda_theta_prio = sparse(lambda_theta_prio);
                
                [J1, d_hatN_dXi, dN_dXi] = obj.compute_jacobian_geometric(Linv, N_control, r, theta, lambda_r_prio, lambda_theta_prio,...
                    phi_r, phi_theta, Nx(i,1), Ny(i,1), Nz(i,1), P_t(1,1), P_t(1,2), P_t(1,3), Phi_t);
                
                jacobian_components.first_term(i).data = J1;
                jacobian_components.second_term(i).data = d_hatN_dXi;
                    [dQ_dk1, dQ_dxi_row, dRho_dk, Lg_rho, Mg_rho, Ng_rho, ... 
                        xi_block] = obj.compute_jacobian_global(d, Linv, N_control, r, theta, lambda_r_prio, lambda_theta_prio,...
                        phi_r, phi_theta, Nx(i,1), Ny(i,1), Nz(i,1), P_t(1,1), P_t(1,2), P_t(1,3), Phi_t, j, control_handles);                
                
                if(mode == 1) % primarily used for global refinement
                    delta_theta1_theta1_ = [delta_theta1_theta1 0];
                    delta_theta1_theta2_ = [delta_theta1_theta2 0];
                    delta_theta2_theta2_ = [delta_theta2_theta2 0];
                    phi_rr = (delta_theta1_theta1_*d) + (Lg_rho*A);
                    phi_r_theta = (delta_theta1_theta2_*d) + (Mg_rho*A); 
                    phi_theta_theta = (delta_theta2_theta2_*d) + (Ng_rho*A);
                    
                    hat_N = [Nx(i,1), Ny(i,1), Nz(i,1), 1];
                    
                    NtN = hat_N.'*hat_N;
                    numer = (phi_rr*NtN*phi_theta_theta.') - (phi_r_theta*NtN*phi_r_theta.');
                    denom_norm = ((Nx(i,1)*N_norm)^2 + (Ny(i,1)*N_norm)^2 + (Nz(i,1)*N_norm)^2);
                    C(i,1) = (numer/denom_norm);
                    
                    %%%% Curvature Jacobian computation %%%
                    xi_block_top = xi_block(1:end-4,:);
                    xi_block_bottom = xi_block(end-3:end,:);
                    d_phi_rrT_d_xi = Lg_rho*xi_block_top;
                    d_phi_rThetaT_d_xi = Mg_rho*xi_block_top;
                    
                    l = size(xi_block_bottom);
                    l = l(2)/4;
                    
                    a1 = xi_block_bottom(1,1:l); a3 = xi_block_bottom(3,1:l);
                    a5 = xi_block_bottom(1,l+1:2*l); a7 = xi_block_bottom(3,l+1:2*l);
                    a9 = xi_block_bottom(1,(2*l)+1:3*l); a11 = xi_block_bottom(3,(2*l)+1:3*l);
                    a13 = xi_block_bottom(1,(3*l)+1:4*l); a15 = xi_block_bottom(3,(3*l)+1:4*l);
                    
                    d_Ng_2_dXi = vertcat(((a1*sT) + (a3*cT)), ((a5*sT) + (a7*cT)), ((a9*sT) + (a11*cT)), ((a13*sT) + (a15*cT)));
                    d_phi_thetaTheta_d_xi = d_phi_rrT_dxi_part1*(kron(Ng_rho.',speye(l))) - d_Ng_2_dXi;
                    
                    d_phi_rTheta_dxi = d_phi_rrT_dxi_part1*(kron(Mg_rho.',speye(l)));
                    
                    dN_dXi_flat = horzcat(d_hatN_dXi(1,:), d_hatN_dXi(2,:),  d_hatN_dXi(3,:),  d_hatN_dXi(4,:));
                    
                    d_N_phiThetaTheta_dXi = ( dN_dXi_flat*( kron(phi_theta_theta.',speye(l)) ) ) + ( hat_N*d_phi_thetaTheta_d_xi ) ;
                    d_phiTrr_N_dXi = ( d_phi_rrT_d_xi*( kron(hat_N.',speye(l)) ) ) + ( phi_rr*d_hatN_dXi ) ;
                    d_N_phiT_rTheta_dXi = ( dN_dXi_flat*( kron(phi_r_theta.',speye(l)) ) ) + ( hat_N*d_phi_rTheta_dxi ) ;
                    d_phi_rTheta_N_dXi = ( d_phi_rThetaT_d_xi*( kron(hat_N.',speye(l)) ) ) + ( phi_r_theta*d_hatN_dXi ) ;
                    
                    d_alphaK_d_xi =  (phi_rr*hat_N.'*d_N_phiThetaTheta_dXi) + ...
                        (d_phiTrr_N_dXi* kron(hat_N*phi_theta_theta.', speye(l)) ) - ...
                        (phi_r_theta*hat_N.'*d_N_phiT_rTheta_dXi) - ...
                        (d_phi_rTheta_N_dXi*( kron(hat_N*phi_r_theta.', speye(l)) ) );
                    
                    dN_dXi_vert = dN_dXi(1:3,:);
                    dN_dXi_horz = horzcat( dN_dXi(1,:) , dN_dXi(2,:) , dN_dXi(3,:) );
                    
                    phiR_cross_phiTheta = hat_N(1,1:3)*N_norm;
                    
                    d_phiR_cross_phiTheta_dXi =  ( dN_dXi_horz * kron(phiR_cross_phiTheta.', speye(l)) ) + ...
                        phiR_cross_phiTheta * dN_dXi_vert;
                    
                    dC_dXi_row = (1/N_norm.^4) * ( ( N_norm.^2 *  d_alphaK_d_xi ) -  ( numer * d_phiR_cross_phiTheta_dXi )  );
                    dC_dXi(i,:) = dC_dXi_row;
                    
                    %%%% Curvature Jacobian computation %%%
                end
                
                jacobian_components.dQ_dk1(i).data = full(dQ_dk1);
                jacobian_components.dRho_dk(i).data = full(dRho_dk);
                
                jacobian_components.dPhi_1(i).data = full(phi_r);
                jacobian_components.dPhi_2(i).data = full(phi_theta);
                jacobian_components.p(i).data = full([r theta]);
                
                if(mode == 0)
                    dQ_dxi(i,:) = dQ_dxi_row;
                else
                    dQ_dxi(i,:) = 0;
                end
                
                Px(i,1) = fp(1,1);
                Py(i,1) = fp(1,2);
                Pz(i,1) = fp(1,3);
                
                if(mode == 0)
                    E = (Px'*xi*Px) + (Py'*xi*Py) + (Pz'*xi*Pz);
                else
                    E = 1;
                end
                
            end
        end
            
            
            %%
            %%%%%%%%%%%%%%%%%%%%
            %  Compute the 3D TPS parameters corresponding to target positions
            %  for the control handles
            %
            % @param: control_handles - the input control handles for the TPS
            % @param: control_handle_target - the target position for the control
            %                                   handles of the TPS
            % @param: N - number of 3D points in data
            %
            % @return: [A, d] - computed parameters of the TPS
            %%%%%%%%%%%%%%%%%%%%
            function [A, d, Linv, xi] = compute_tps(obj, control_handles, control_handle_target, N)
                K = zeros(N, N, 'double');
                P = ones (N, 4, 'double');
                O = zeros(4, 4, 'double');
                Y = zeros((N+4), 4, 'double');
                
                for i = 1:N
                    for j = 1:N
                        if(i~=j)
                            x_n = control_handles(i,1) - control_handles(j,1);
                            y_n = control_handles(i,2) - control_handles(j,2);
                            z_n = control_handles(i,3) - control_handles(j,3);
                            K(i,j) = obj.U(x_n, y_n, z_n);
                        end
                    end
                    P(i,1) = control_handles(i,1);
                    P(i,2) = control_handles(i,2);
                    P(i,3) = control_handles(i,3);
                    P(i,4) = 1.0;
                    Y(i,1) = control_handle_target(i,1);
                    Y(i,2) = control_handle_target(i,2);
                    Y(i,3) = control_handle_target(i,3);
                    Y(i,4) = 1.0;
                end
                
                Kinv = pinv(K);
                xi = Kinv - ( Kinv *  P * pinv(P' * Kinv * P) * P' * Kinv );
                
                KP = horzcat(K,P);
                PTO = horzcat(P',O);
                L = vertcat(KP,PTO);
                
                Linv = pinv(L);
                W = Linv*Y;
                A = W(1:N,:);
                d = W(N+1:end,:);
            end
            
            %% compute Jacobian geometric
            function [J1, dN_dY, dN_dY_1] = compute_jacobian_geometric(obj, Linv, N_control, theta1, theta2,...
                    lambda_r_prio, lambda_theta_prio, phi_r, phi_theta,...
                    Nx, Ny, Nz, Px, Py, Pz, phi)
                theta1 = deg2rad(theta1);            
                theta2 = deg2rad(theta2);            
                ones_block = speye(N_control, N_control);
                ones_block_thrice = speye(3*N_control, 3*N_control);
                zeros_block = zeros(N_control, (3 * N_control),'like',ones_block);
                zeros_row = zeros(1, (3 * N_control),'like',ones_block);
                block_matrix_top = horzcat(ones_block, zeros_block, ones_block, zeros_block, ones_block, zeros_block);
                block_matrix_bottom = zeros(4,(12*N_control),'like',ones_block);
                block_matrix_total = vertcat(block_matrix_top, block_matrix_bottom);
                
                bar_xi_lambda = Linv([1:N_control],[1:N_control]);
                Linv_bottom_left = Linv([N_control+1 : N_control + 4],[1 : N_control]);
                angular_component1 = [1 0 0];
                angular_component2 = [0 1 0];
                
                
                delta_a_delta_Y = Linv_bottom_left * block_matrix_top;
                
                delA1 = delta_a_delta_Y(1,:);
                delA2 = delta_a_delta_Y(2,:);
                delA3 = delta_a_delta_Y(3,:);
                
                delA_theta = vertcat(delA1, delA2, delA3);
                
                
                del_r_del_Y = angular_component1*delA_theta + (lambda_r_prio * bar_xi_lambda * block_matrix_top);
                
                
                del_theta_del_Y_1 = angular_component2*delA_theta;
                del_theta_del_Y_2 = lambda_theta_prio * bar_xi_lambda * block_matrix_top;
                
                del_theta_del_Y = del_theta_del_Y_1 + del_theta_del_Y_2;
                
                phi_r_skew = [0 -phi_r(3) phi_r(2) ; phi_r(3) 0 -phi_r(1) ; -phi_r(2) phi_r(1) 0 ];
                
                
                del_r_del_Y_x = del_r_del_Y(1:(3*N_control));
                del_r_del_Y_y = del_r_del_Y((3*N_control+1):(6*N_control));
                del_r_del_Y_z = del_r_del_Y((6*N_control+1):(9*N_control));
                
                dNdY_block1_1 = horzcat(zeros_row, -del_r_del_Y_z, del_r_del_Y_y);
                dNdY_block1_2 = horzcat(del_r_del_Y_z, zeros_row, -del_r_del_Y_x);
                dNdY_block1_3 = horzcat(-del_r_del_Y_y, del_r_del_Y_x, zeros_row);
                
                dNdY_block1 = vertcat(dNdY_block1_1, dNdY_block1_2, dNdY_block1_3);
                dNdY_block2 = kron(phi_theta(1,1:3).', ones_block_thrice);
                
                del_theta_del_Y_x = del_theta_del_Y(1:(3*N_control));
                del_theta_del_Y_y = del_theta_del_Y((3*N_control+1):(6*N_control));
                del_theta_del_Y_z = del_theta_del_Y((6*N_control+1):(9*N_control));
                
                dNdY_block3 = vertcat(del_theta_del_Y_x, del_theta_del_Y_y, del_theta_del_Y_z);
                
                dN_tilde_dY = (dNdY_block1 * dNdY_block2) + (phi_r_skew * dNdY_block3);
                
                norm_N = sqrt(Nx.^2 + Ny.^2 + Nz.^2);
                norm_N3 = norm_N.^3;
                
                dN_dY_1 = vertcat(dN_tilde_dY, zeros_row);
                dN_dY_1 = dN_dY_1/norm_N;
                
                dN_dY_2_1 = (Nx*dN_tilde_dY(1,:) + Ny*dN_tilde_dY(2,:) + Nz*dN_tilde_dY(3,:));
                
                dN_dY_2 = vertcat((Nx*dN_dY_2_1), (Ny*dN_dY_2_1), (Nz*dN_dY_2_1), (dN_dY_2_1));
                
                dN_dY = dN_dY_1 - (dN_dY_2/norm_N3);
                
                P = [Px Py Pz 1];
                
                du_dY_1 = horzcat(phi,P);
                
                du_dY = du_dY_1 * Linv * block_matrix_total;
                
                N_j = [Nx;Ny;Nz;1];
                
                J1 = du_dY* kron(N_j,ones_block_thrice);
            end
            
            %%  compute jacobian global refinement
            function [dQ_dk1, dQ_dxi, dRho_dk, Lg_rho, Mg_rho, Ng_rho, ... 
                    xi_block] = compute_jacobian_global(obj, D, Linv, N_control, r, theta,...
                    lambda_r_prio, lambda_theta_prio, phi_r, phi_theta,...
                    Nx, Ny, Nz, Px, Py, Pz, phi, j, control_handles)
                
                x = r;
                y = theta;
                [X_, Y_, Z_] = obj.cartesian_coordinates(x, y);
                
                dx_dk = zeros(1, (2*N_control));
                dy_dk = zeros(1, (2*N_control));
                
                ones_block = speye(N_control, N_control);
                zeros_block = zeros(N_control, (3 * N_control),'like',ones_block);
                zeros_row = zeros(1, (2 * N_control),'like',ones_block);
                block_matrix_top = horzcat(ones_block, zeros_block, ones_block, zeros_block, ones_block, zeros_block);
                block_matrix_bottom = zeros(4,(12*N_control));
                block_matrix_total = vertcat(block_matrix_top, block_matrix_bottom);
                block_matrix_total = sparse(block_matrix_total);
                P = [Px Py Pz 1];
                
                dQ_dxi_1 = horzcat(phi, P);
                
                Linv = sparse(Linv);
                
                xi_block =  Linv * block_matrix_total;
                dQ_dxi = dQ_dxi_1 * xi_block;
                
                dRho_dk = zeros(1, (2*N_control*N_control),'like',ones_block);
                
                Lg_rho = zeros(1, N_control);
                Mg_rho = zeros(1, N_control);
                Ng_rho = zeros(1, N_control);
                
                
                P_canonical = [X_ Y_ Z_];
                
                for i = 1:N_control % -> possible issue with indices
                    dRho_dk(1, ((i-1)*N_control) + j) =   (x - control_handles(j,1)) /(phi(1,i).^2);
                    dRho_dk(1, ((i-1)*N_control) + (N_control + j)) = (y - control_handles(j,1)) /(phi(1,i).^2);
                    
                    diff_C = (control_handles(i,:) - P_canonical);
                    Lg_rho(1,i) = ( phi(1,i).^2  -  diff_C(1,1).^2 ) / phi(1,i).^2;
                    Mg_rho(1,i) = ( phi(1,i).^2  -  ( diff_C(1,1) + diff_C(1,2)).^2 ) / phi(1,i).^2;
                    Ng_rho(1,i) = ( phi(1,i).^2  -  diff_C(1,2).^2 ) / phi(1,i).^2;
                    
                end
                Lg_rho = sparse(Lg_rho);
                Mg_rho = sparse(Mg_rho);
                Ng_rho = sparse(Ng_rho);
                
                
                dx_dk(1,j) = 1;
                dy_dk(1, N_control + j) = 1;
                
                dRho_dk = sparse(dRho_dk);
                dx_dk = sparse(dx_dk);
                dy_dk = sparse(dy_dk);
                dQ_dk1 = horzcat(dRho_dk, dx_dk, dy_dk, zeros_row, zeros_row);
                
            end
            
            %%
            %%%%%%%%%%%%%%%%%%%%
            %  Update the TPS based on updated positions of the control handle
            %
            % @param: control_handle_target - the updated target position for the control
            %                                   handles of the TPS
            % @param: Linv - inverse of design matrix, pre-computed
            % @param: xi - bending energy matrix
            %
            % @return: A, d - updated computed parameters of the TPS
            %                  E, J - bending energy evaluated at the control
            %                  handles and their Jacobian
            %%%%%%%%%%%%%%%%%%%%
            function [A, d, E,J] = update_tps(obj, control_handle_target, Linv, xi, N, quickReturn)
                Y = zeros((N+4), 4, 'double');
                
                for i = 1:N
                    Y(i,1) = control_handle_target(i,1);
                    Y(i,2) = control_handle_target(i,2);
                    Y(i,3) = control_handle_target(i,3);
                    Y(i,4) = 1.0;
                end
                
                W = Linv*Y;
                A = W(1:N,:);
                d = W(N+1:end,:);

                if exist('quickReturn', 'var')
                    E = []; J = [];
                    return;
                end
                
                E = Y(1:N,1).'*xi*Y(1:N,1) + Y(1:N,2).'*xi*Y(1:N,2) + Y(1:N,3).'*xi*Y(1:N,3);
                J = zeros(1,(3*N));
                
                dYxT_dY = zeros(1,(3*N*N));
                dYyT_dY = zeros(1,(3*N*N));
                dYzT_dY = zeros(1,(3*N*N));

                ones_block = speye(N,N);
                zeros_block = zeros(N,N);
                ones_block_thrice = speye(3*N,3*N);

                dYx_dY = horzcat(ones_block, zeros_block, zeros_block);
                dYy_dY = horzcat(zeros_block, ones_block, zeros_block);
                dYz_dY = horzcat(zeros_block, zeros_block, ones_block);

                step = 0;
                for i = 1:N
                    dYxT_dY(1,i + (step*N*3)) = 1;
                    dYyT_dY(1,i + N + (step*N*3)) = 1;
                    dYzT_dY(1,i + N + N + (step*N*3)) = 1;
                    step = step + 1;
                end

                J(1,:) = ( dYxT_dY * kron((xi*Y(1:N,1)), ones_block_thrice) + Y(1:N,1).'*xi*dYx_dY + ...
                            dYyT_dY * kron((xi*Y(1:N,2)), ones_block_thrice) + Y(1:N,2).'*xi*dYy_dY + ...
                            dYzT_dY * kron((xi*Y(1:N,3)), ones_block_thrice) + Y(1:N,3).'*xi*dYz_dY );
                
            end
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Generate a sphere of uniformly spaced vertices which could be used as control handles
            %
            % @param: num_points - the number of points per axis of the sphere will be approximately equal to sqrt{num_points}
            % @return: [control_handles] - 3D coordinates of the vertices of the sphere
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [control_handles] = generate_control_handles(obj, num_points)
                num_pts_per_axis = round(sqrt(num_points));
                
                theta = linspace(0, 2*pi, num_pts_per_axis);
                X = ones(num_pts_per_axis,1); Y = ones(num_pts_per_axis,1); Z = ones(num_pts_per_axis,1);
                num_point = 1;
                for index1 = 1:num_pts_per_axis
                    for index2 = 1:num_pts_per_axis
                        [X(num_point,1), Y(num_point,1), Z(num_point,1)] = obj.cartesian_coordinates(theta(index1), theta(index2));
                        num_point = num_point + 1;
                    end
                end
                control_handles = horzcat(X, Y, Z);
            end
            
          %%
            %%%%%%%%%%%%%%%%%%%%
            %  U value for the TPS
            %  CAUTION: This is just a L2 norm function, and not the complete U value
            %  The call to this function needs to handle generating the U value
            %
            % @param: x,y, z - cartesian coordinates of the point
            %
            % @return: Ur - computed U-value for this point
            %%%%%%%%%%%%%%%%%%%%
            function Ur = U(obj, x,y,z)
                Ur = sqrt((x*x) + (y*y) + (z*z));
                if(abs(Ur) < 0.0000001)
                    Ur = 0.0001;
                end
            end
            
          %%
            %%%%%%%%%%%%%%%%%%%%
            %  Find the normal to the parameterized surface
            %
            % @param: A - A matrix corresponding to the TPS
            % @param: D - D matrix corresponding to the TPS
            % @param: control_handles - the input control handles for the TPS
            % @param: r - the r coordinate of the point for which the normal needs to
            %             be computed
            % @param: theta - the theta coordinate of the point for which the
            %              normal needs to be computed
            % @param: N - number of 3D points in data
            %
            % @return: [Nx, Ny, Nz] - the computed normals
            %%%%%%%%%%%%%%%%%%%%
            function [Nx, Ny, Nz,lambda_r_prio, lambda_theta_prio, phi_r, phi_theta, norm, cT, sT] = get_normal(obj, A, D, control_handles, phi, theta, X, Y, Z, N)
                
                Ut_phi = zeros(1,N,'like',speye(2));
                Ut_theta = zeros(1,N,'like',speye(2));
                
                D_coeff_phi = zeros(1,4,'like',speye(2));
                D_coeff_theta = zeros(1,4,'like',speye(2));
                
                M_phi     = zeros(1,3,'like',speye(2));
                M_theta   = zeros(1,3,'like',speye(2));
                
                
                D_coeff_phi(1,1) = 1;  % Changes needed in Line 70-71 as well
                D_coeff_phi(1,2) = 0;
                D_coeff_phi(1,3) = 0;
                D_coeff_phi(1,4) = 0;
                
                D_coeff_theta(1,1) = 0;
                D_coeff_theta(1,2) = 1;
                D_coeff_theta(1,3) = 0;
                D_coeff_theta(1,4) = 0;
                
                
                for index = 1:N
                    Xc = control_handles(index,1);
                    Yc = control_handles(index,2);
                    Zc = control_handles(index,3);
                    
                    Ut = obj.U((X - Xc),(Y - Yc),(Z - Zc));
                    if(abs(Ut) > 0.000001)
                        Ut_phi(1,index) = (1/Ut)*(X - Xc);
                        Ut_theta(1,index) = (1/Ut)*(Y - Yc);
                    else
                        Ut_phi(1,index) = 0.0;
                        Ut_theta(1,index) = 0.0;
                    end
                end
                
                
                Phi_phi = (D_coeff_phi*D) + (Ut_phi*A);
                Phi_theta = (D_coeff_theta*D) + (Ut_theta*A);
                
                M_phi(1,1) = Phi_phi(1,1); M_phi(1,2) = Phi_phi(1,2); M_phi(1,3) = Phi_phi(1,3);
                M_theta(1,1) = Phi_theta(1,1); M_theta(1,2) = Phi_theta(1,2); M_theta(1,3) = Phi_theta(1,3);
                
                N = cross(M_phi,M_theta);
                normN = sqrt((N(1,1)*N(1,1)) + (N(1,2)*N(1,2)) + (N(1,3)*N(1,3)));
                
                N = N / normN;
                
                Nx = N(1,1);
                Ny = N(1,2);
                Nz = N(1,3);
                
                %%%Populating the output variables in sync with WTT format
                lambda_r_prio = Ut_phi;
                lambda_theta_prio = Ut_theta;
                phi_r = M_phi;
                phi_theta = M_theta;
                norm = normN;
                cT = cos(theta); sT = sin(theta);
                
            end

        end
    end
