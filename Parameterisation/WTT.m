
classdef WTT %< handle
    
    methods
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Interpolate the 3D points from undeformed, right circular cylinder     %
        %  based on previously obtained TPS values                                %
        %                                                                         %
        % @param: A, d - TPS parameter                                            %
        % @param: interpolation_candidate_vertices - vertices on right circular   %
        %         cylinder that needs to be interpolated, could be of higher resolution                          %
        % @param: control_handles - undeformed/original control handles defined   %
        %                           for the TPS                                   %
        % @param: N_test - number of vertices in the 'interpolation_candidate_vertices' (which will be     %
        %                     interpolated)                                       %
        % @param: N_control - number of 3D points in 'control_handles'        %
        %                                                                         % 
        % @param: jacobian_mode - 0 -> initial refinement; 1 -> global
        % refinement
        % @return: [Px, Py, Pz, Nx, Ny, Nz] - interpolated points & normals       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Px, Py, Pz, Nx, Ny, Nz, jacobian_components] = interpolate_tps_simple(obj, A, d, interpolation_candidate_vertices, ...
                                                                    control_handles, N_test, N_control, xi, Linv, jacobian_mode)



            Px =  zeros(N_test, 1, 'double');
            Py =  zeros(N_test, 1, 'double');
            Pz =  zeros(N_test, 1, 'double');
            Nx = zeros(N_test, 1, 'double');
            Ny = zeros(N_test, 1, 'double');
            Nz = zeros(N_test, 1, 'double');
%             E =   zeros(N_test, 1, 'double');

            index = 1;

            for i = 1:N_test
               P_t = ones(1, 4, 'double');
               P_t(1,1) = interpolation_candidate_vertices(i,1);
               P_t(1,2) = interpolation_candidate_vertices(i,2);
               P_t(1,3) = interpolation_candidate_vertices(i,3);
               Phi_t = zeros(1, N_control, 'double');
               
               for j = 1:N_control
                   Phi_t(1,j) = obj.U(  control_handles(j,1) -   P_t(1,1),     ...
                                        control_handles(j,2) -   P_t(1,2),     ...
                                        control_handles(j,3) -   P_t(1,3));
               end
               
               fp = (P_t*d) + (Phi_t*A);       

               [r, theta] = obj.flattened_coordinates(P_t(1,1), P_t(1,2), P_t(1,3));
               [Nx(i,1), Ny(i,1), Nz(i,1), lambda_r_prio, lambda_theta_prio, phi_r, phi_theta] = obj.get_normal(A, d, control_handles, r, theta, N_control);
%                [E(i,1), phi_rr, phi_tt, phi_rt] = obj.get_bending_energy(A, d, control_handles, r, theta, N_control);
               
               if(jacobian_mode == 0) % computing jacobian stuff for init. parameterisation
                    [J1, dN_dY] = obj.compute_jacobian_geometric(Linv, N_control, theta, lambda_r_prio, lambda_theta_prio,...
                        phi_r, phi_theta, Nx(i,1), Ny(i,1), Nz(i,1), P_t(1,1), P_t(1,2), P_t(1,3), Phi_t);
               
                    jacobian_components.first_term(i).data = J1;
                    jacobian_components.second_term(i).data = dN_dY;
               else % computing jacobian stuff for global refinement
                   
%                    fprintf('FATAL ERROR: worng option for type of Jacobian required! Code will crash soon.');
                    jacobian_components = 0;
               end
               
%                E(i,1) = 0;
%                [E_] = obj.get_bending_energy_3D(A, control_handles, P_t(1,1), P_t(1,2), P_t(1,3), ... 
%                                                                 Phi_t, N_control);
%                E(index,1) = E_(1,1); index = index + 1;
%                E(index,1) = E_(1,2); index = index + 1;
%                E(index,1) = E_(1,3); index = index + 1;
               Px(i,1) = fp(1,1);
               Py(i,1) = fp(1,2);
               Pz(i,1) = fp(1,3);
            end
            
%             E = (Px'*xi*Px) + (Py'*xi*Py) + (Pz'*xi*Pz);
        end
        
                %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Interpolate the 3D points from undeformed, right circular cylinder     %
        %  based on previously obtained TPS values                                %
        %                                                                         %
        % @param: A, d - TPS parameter                                            %
        % @param: interpolation_candidate_vertices - vertices on right circular   %
        %         cylinder that needs to be interpolated                          %
        % @param: control_handles - undeformed/original control handles defined   %
        %                           for the TPS                                   %
        % @param: N_test - number of vertices in the test data (which will be     %
        %                     interpolated)                                       %
        % @param: N_control - number of 3D points in data                         %
        %                                                                         % 
        % @param: jacobian_mode - 0 -> isometric OR reprojective error; 1 -> smoothing error with high-res mesh
        % refinement
        % @return: [Px, Py, Pz, Nx, Ny, Nz] - interpolated points & normals       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Px, Py, Pz, Nx, Ny, Nz, E, jacobian_components, dQ_dxi, C, dC_dXi] = interpolate_tps(obj, A, d, template, ...
                                                                    source_points, N_test, N_control, xi, Linv, jacobian_mode)

            % jacobian_mode: -1 : just points, nothing else
            % jacobian_mode: -2: just points & normals, nothing else
            if(jacobian_mode<=-1)
                Px =  zeros(N_test, 1);
                Py =  zeros(N_test, 1);
                Pz =  zeros(N_test, 1);
                Nx = zeros(N_test, 1);
                Ny = zeros(N_test, 1);
                Nz = zeros(N_test, 1);
                for i = 1:N_test
                    P_t = ones(1, 4, 'double');
                    P_t(1,1) = template(i,1);
                    P_t(1,2) = template(i,2);
                    P_t(1,3) = template(i,3);
                    Phi_t = zeros(1, N_control);

                    for j = 1:N_control
                        Phi_t(1,j) = obj.U(  source_points(j,1) -  P_t(1,1),  ...
                            source_points(j,2) -  P_t(1,2), ...
                            source_points(j,3) -  P_t(1,3)  );
                    end
                    fp = (P_t*d) + (Phi_t*A);
                    if(jacobian_mode == -2)
                        [r, theta] = obj.flattened_coordinates(P_t(1,1), P_t(1,2), P_t(1,3));
    
                        [Nx(i,1), Ny(i,1), Nz(i,1), ~, ~, ~, ~, ~, ~, ~] = obj.get_normal(A, d, source_points, r, theta, N_control);
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

            Px =  zeros(N_test, 1, 'like', speye(3));
            Py =  zeros(N_test, 1, 'like', speye(3));
            Pz =  zeros(N_test, 1, 'like', speye(3));
            Nx = zeros(N_test, 1, 'like', speye(3));
            Ny = zeros(N_test, 1, 'like', speye(3));
            Nz = zeros(N_test, 1, 'like', speye(3));
            E =   zeros(N_test, 1, 'like', speye(3));
            C = zeros(N_test,1, 'like', speye(3));
%             jacobian_components = 0; %%%-> remove this!
            
            dQ_dxi = zeros(N_test, (12*N_test), 'like', speye(3));
            
            dC_dXi = zeros(N_test, (3*N_control), 'like', speye(3));

            index = 1;
            dYT_dxi = zeros(4,(3*N_control), 'like', speye(3));
            offset = 0;
            for loop = 1:(N_control)
                dYT_dxi(1,loop+(3*offset*N_control)) = 1.0;
                dYT_dxi(2,loop+(3*offset*N_control) + N_control) = 1.0;
                dYT_dxi(3,loop+(3*offset*N_control) + N_control + N_control) = 1.0;
                offset = offset + 1;
            end
            
            I_3l = speye(3*N_control);
            d_phi_rrT_dxi_part1 = dYT_dxi*(kron(xi, I_3l));

            for i = 1:N_test
               P_t = ones(1, 4, 'double');
               
% % % %                r = template(i,2);
% % % %                theta = template(i,1);
% % % %                theta = rad2deg(theta);
% % % %                
% % % %                
% % % %                P_t(1,1) = sin(theta);
% % % %                P_t(1,2) = r;
% % % %                P_t(1,3) = cos(theta);

               P_t(1,1) = template(i,1);
               P_t(1,2) = template(i,2);
               P_t(1,3) = template(i,3);
               [r, theta] = obj.flattened_coordinates(P_t(1,1), P_t(1,2), P_t(1,3));
               
% % % %                fprintf('%0.4f %0.4f %0.4f\n',P_t(1,1),P_t(1,2),P_t(1,3));
% % % %                fprintf('%0.4f \n',theta);

               Phi_t = zeros(1, N_control, 'like', speye(3));
              
               
               for j = 1:N_control
                   Phi_t(1,j) = obj.U( source_points(j,1) -   P_t(1,1),     ...
                                                source_points(j,2) -   P_t(1,2),     ...
                                                source_points(j,3) -   P_t(1,3));
               end
               
               Phi_t = sparse(Phi_t);
               A = sparse(A);
               
               fp = (P_t*d) + (Phi_t*A);       
% % % %                fprintf('%0.4f %0.4f %0.4f\n',fp(1,1),fp(1,2),fp(1,3));

% % % %                [r, theta] = obj.cylindrical_coordinates_Yaligned_forward(P_t(1,1), P_t(1,2), P_t(1,3));
               [Nx(i,1), Ny(i,1), Nz(i,1), lambda_r_prio, lambda_theta_prio, phi_r, phi_theta, N_norm, cT, sT] = obj.get_normal(A, d, source_points, r, theta, N_control);
% % % %                [E(i,1), phi_rr, phi_tt, phi_rt] = obj.get_bending_energy(A, d, control_handles, r, theta, N_control);

                    [J1, d_hatN_dXi, dN_dXi] = obj.compute_jacobian_geometric(Linv, N_control, theta, lambda_r_prio, lambda_theta_prio,...
                        phi_r, phi_theta, Nx(i,1), Ny(i,1), Nz(i,1), P_t(1,1), P_t(1,2), P_t(1,3), Phi_t);
                
                jacobian_components.first_term(i).data = J1;
                jacobian_components.second_term(i).data = d_hatN_dXi;
               
                [dQ_dk1, dQ_dxi_row, dRho_dk, dSin_dk, dR_dk, dCos_dk, Lg_rho, Mg_rho, Ng_rho, Ng_2, xi_block] = obj.compute_jacobian_global(d, Linv, N_control, r, theta, lambda_r_prio, lambda_theta_prio,...
                    phi_r, phi_theta, Nx(i,1), Ny(i,1), Nz(i,1), P_t(1,1), P_t(1,2), P_t(1,3), Phi_t, j, source_points);
                
                phi_rr = Lg_rho*A;
                phi_r_theta = Mg_rho*A;
                Ng_1_1 = Ng_rho*A;
                phi_theta_theta = Ng_1_1 - Ng_2;
                
                hat_N = [Nx(i,1), Ny(i,1), Nz(i,1), 1];
                
                NtN = hat_N.'*hat_N;
                numer = (phi_rr*NtN*phi_theta_theta.') - (phi_r_theta*NtN*phi_r_theta.');
                denom_norm = ((Nx(i,1)*N_norm)^2 + (Ny(i,1)*N_norm)^2 + (Nz(i,1)*N_norm)^2);
                C(i,1) = (numer/denom_norm);
                
                %%%% test code
% % % %                 NtN = hat_N*hat_N.';
% % % %                 numer = (phi_rr*NtN*phi_theta_theta.') - (phi_r_theta*NtN*phi_r_theta.');
% % % %                 denom_norm = ((Nx(i,1)*N_norm)^2 + (Ny(i,1)*N_norm)^2 + (Nz(i,1)*N_norm)^2);
% % % %                 C_test1 = (numer/denom_norm);
% % % %                 
% % % %                 L = dot(phi_rr, hat_N);
% % % %                 M = dot(phi_r_theta, hat_N);
% % % %                 N = dot(phi_theta_theta, hat_N); % problem with either L, M or N
% % % %                 
% % % %                 E_ = (phi_r(1,1).^2 + phi_r(1,2).^2 + phi_r(1,3).^2);
% % % %                 F_  =  (phi_r(1,1)*phi_theta(1,1)) + (phi_r(1,2)*phi_theta(1,2)) + (phi_r(1,3)*phi_theta(1,3));
% % % %                 G_ = (phi_theta(1,1).^2 + phi_theta(1,2).^2 + phi_theta(1,3).^2);
% % % %                 
% % % %                 C_test2 = ((L*N) - (M*M)) / ((E_*G_) - (F_*F_));
% % % %                 
% % % %                 fprintf('%d %d %d\n',C(i,1), C_test1, C_test2);
                %%% test code ends
                
                
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
                
%                 d_Ng_2_dXi = horzcat(((a1*sT) + (a3*cT)), ((a5*sT) + (a7*cT)), ((a9*sT) + (a11*cT)), ((a13*sT) + (a15*cT)));                
%                 d_phi_thetaThetaT_d_xi =  (Ng_rho*xi_block_top) - d_Ng_2_dXi;

                 d_Ng_2_dXi = vertcat(((a1*sT) + (a3*cT)), ((a5*sT) + (a7*cT)), ((a9*sT) + (a11*cT)), ((a13*sT) + (a15*cT)));
                d_phi_thetaTheta_d_xi = d_phi_rrT_dxi_part1*(kron(Ng_rho.',speye(l))) - d_Ng_2_dXi;
                
%                 d_phi_rr_dxi = d_phi_rrT_dxi_part1*(kron(Lg_rho.',eye(l)));
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
% % % %                 
% % % %                 fprintf('Curvature:\t%0.4f\t%0.4f\n',C(i,1),(numer/denom_norm));

                jacobian_components.dQ_dk1(i).data = dQ_dk1;
                jacobian_components.dRho_dk(i).data = dRho_dk;
                jacobian_components.dSin_dk(i).data = dSin_dk;
                jacobian_components.dR_dk(i).data = dR_dk;
                jacobian_components.dCos_dk(i).data = dCos_dk;
                
                jacobian_components.dPhi_1(i).data = phi_r;
                jacobian_components.dPhi_2(i).data = phi_theta;
                jacobian_components.p(i).data = [r theta];
                
                if(jacobian_mode == 0)
                    dQ_dxi(i,:) = dQ_dxi_row; 
                else
                    dQ_dxi(i,:) = 0;
                end

               Px(i,1) = fp(1,1); % -> uncomment
               Py(i,1) = fp(1,2);
               Pz(i,1) = fp(1,3);               
            end
            
            if(jacobian_mode == 0)
                E = (Px'*xi*Px) + (Py'*xi*Py) + (Pz'*xi*Pz);
            else
                E = 1;
            end
            
        end
        
%% compute jacobian global refinement
function [dQ_dk1, dQ_dxi, dRho_dk, dSin_dk, dR_dk, dCos_dk, Lg_rho, Mg_rho, Ng_rho, Ng_2, xi_block] = compute_jacobian_global(obj, D, Linv, N_control, r, theta,... 
                                                                        lambda_r_prio, lambda_theta_prio, phi_r, phi_theta,... 
                                                                        Nx, Ny, Nz, Px, Py, Pz, phi, j, control_handles)
                                                                    
        theta = deg2rad(theta);

        ones_block = speye(N_control, N_control);
        zeros_block = zeros(N_control, (3 * N_control), 'like', speye(3));
        zeros_row = zeros(1, (2 * N_control), 'like', speye(3));
        block_matrix_top = horzcat(ones_block, zeros_block, ones_block, zeros_block, ones_block, zeros_block);
        block_matrix_bottom = zeros(4,(12*N_control), 'like', speye(3));
        block_matrix_total = vertcat(block_matrix_top, block_matrix_bottom);
        P = [Px Py Pz 1];
        
        a1 = D(1,1); a5 = D(1,2); a9  = D(1,3); a13  = D(1,4);
        a3 = D(3,1); a7 = D(3,2); a11 = D(3,3); a15  = D(2,4);
        
        dQ_dxi_1 = horzcat(phi, P);
        
        xi_block =  Linv * block_matrix_total;
        xi_block = sparse(xi_block);
        dQ_dxi = dQ_dxi_1 * xi_block;
        
        dRho_dk = zeros(1, (2*N_control*N_control), 'like', speye(3));
        dSin_dk = zeros(1, (2*N_control), 'like', speye(3));
        dR_dk = zeros(1, (2*N_control), 'like', speye(3));
        dCos_dk = zeros(1, (2*N_control), 'like', speye(3));
        
        Lg_rho = zeros(1, N_control, 'like', speye(3));
        Mg_rho = zeros(1, N_control, 'like', speye(3));
        Ng_rho = zeros(1, N_control, 'like', speye(3));
        Ng_2 = zeros(1, 4, 'like', speye(3));
        
%         fprintf('%d %d\n',r,theta);
        
        for i = 1:N_control % -> possible issue with indices
%             fprintf('\n\n[WTT.m:compute_jacobian_global] -> BIG PROBLEM with dRho_dk (j is a constant index), FIX needed!!!\n\n')
            dRho_dk(1, ((i-1)*N_control) + j) = -((r - control_handles(j,2))/phi(1,i));
            dRho_dk(1, ((i-1)*N_control) + (N_control + j)) = -( ((control_handles(j,3) * sin(theta)) - (control_handles(j,1) * cos(theta)))/phi(1,i) );
            
            Cx = control_handles(i,1);
            Cy = control_handles(i,2);
            Cz = control_handles(i,3);
            
            Lg_rho(1,i) = ( phi(1,i).^2 - (r - Cy).^2 )/( phi(1,i).^3 );
            Mg_rho(1,i) = ( - ( (r - Cy)*( (Cz * sin(theta)) - (Cx * cos(theta)) ) ) )/( phi(1,i).^3 );
            Ng_rho(1,i) =  ( phi(1,i).^2*( (Cz*cos(theta)) + (Cx*sin(theta)) )... 
                                            - ( (Cz * sin(theta)) - (Cx * cos(theta)) ).^2 )/ (phi(1,i).^3);            
        end
        
        Lg_rho = sparse(Lg_rho);
        Mg_rho = sparse(Mg_rho);
        Ng_rho = sparse(Ng_rho);
        
        Ng_2(1,1) = (a1*sin(theta)) + (a3*cos(theta));
        Ng_2(1,2) = (a5*sin(theta)) + (a7*cos(theta));
        Ng_2(1,3) = (a9*sin(theta)) + (a11*cos(theta));
        Ng_2(1,4) = (a13*sin(theta)) + (a15*cos(theta));

        dSin_dk(1, N_control + j) = cos(theta);
        dR_dk(1, j) = 1;
        dCos_dk(1, N_control + j) = -sin(theta);
        
        dQ_dk1 = horzcat(dRho_dk, dSin_dk, dR_dk, dCos_dk, zeros_row);
        dQ_dk1 = sparse(dQ_dk1);
        
end
        
        
%% compute Jacobian geometric
    function [J1, dN_dY, dN_dY_1] = compute_jacobian_geometric(obj, Linv, N_control, theta,... 
                                                                        lambda_r_prio, lambda_theta_prio, phi_r, phi_theta,... 
                                                                                                                        Nx, Ny, Nz, Px, Py, Pz, phi)
        theta = deg2rad(theta);                                                                                                                    
        ones_block = speye(N_control, N_control);
        ones_block_thrice = speye(3*N_control, 3*N_control);
        zeros_block = zeros(N_control, (3 * N_control), 'like', speye(3));
        zeros_row = zeros(1, (3 * N_control), 'like', speye(3));
        block_matrix_top = horzcat(ones_block, zeros_block, ones_block, zeros_block, ones_block, zeros_block);
        block_matrix_bottom = zeros(4,(12*N_control), 'like', speye(3));
        block_matrix_total = vertcat(block_matrix_top, block_matrix_bottom);
        
        block_matrix_top = sparse(block_matrix_top);
        block_matrix_total = sparse(block_matrix_total);
        
        bar_xi_lambda = Linv([1:N_control],[1:N_control]);
        Linv_bottom_left = Linv([N_control+1 : N_control + 4],[1 : N_control]);
        angular_component = [ cos(theta)  -sin(theta) ];
        
        
        delta_a_delta_Y = Linv_bottom_left * block_matrix_top;
        
        delA1 = delta_a_delta_Y(1,:);
        delA2 = delta_a_delta_Y(2,:);
        delA3 = delta_a_delta_Y(3,:);
        
%         a1 = delta_a_delta_Y(1,1:(3*N_control)); a5 = delta_a_delta_Y(1,(3*N_control+1):(6*N_control)); a9 = delta_a_delta_Y(1,(6*N_control+1):(9*N_control)); a13 = delta_a_delta_Y(1,(9*N_control+1):(12*N_control));
%         a3 = delta_a_delta_Y(3,1:(3*N_control)); a7 = delta_a_delta_Y(3,(3*N_control+1):(6*N_control)); a11 = delta_a_delta_Y(3,(6*N_control+1):(9*N_control)); a15 = delta_a_delta_Y(3,(9*N_control+1):(12*N_control));
        
        delA_theta = vertcat(delA1, delA3);
        
%         angular_component_augmented = kron(angular_component, eye((3*N_control),(3*N_control)));
        
        
        del_r_del_Y = delA2 + (lambda_r_prio * bar_xi_lambda * block_matrix_top);
        
        
%         del_theta_del_Y_1 = delA_theta * angular_component_augmented;
        del_theta_del_Y_1 = angular_component*delA_theta;
        del_theta_del_Y_2 = lambda_theta_prio * bar_xi_lambda * block_matrix_top;
        
        del_theta_del_Y = del_theta_del_Y_1 + del_theta_del_Y_2;
        
        phi_r_skew = [0 -phi_r(3) phi_r(2) ; phi_r(3) 0 -phi_r(1) ; -phi_r(2) phi_r(1) 0 ];
%         phi_theta_skew = [0 -phi_theta(3) phi_theta(2) ; phi_theta(3) 0 -phi_theta(1) ; -phi_theta(2) phi_theta(1) 0 ];


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

%%
        function [E] = get_bending_energy_3D(obj, A, control_handles, X, Y, Z, U_tp, N)           
            d2Utp_dX2 = zeros(1,N);
            d2Utp_dY2 = zeros(1,N);
            d2Utp_dZ2 = zeros(1,N);
            
            d2Utp_dXY = zeros(1,N);
            d2Utp_dXZ = zeros(1,N);
            d2Utp_dYZ = zeros(1,N);
            
            for index = 1:N
                Ut = U_tp(1,N);
                Xc = control_handles(index,1);
                Yc = control_handles(index,2);
                Zc = control_handles(index,3);

                d2Utp_dX2(1,index) = (Ut - ((X - Xc)/(2*Ut)) )/(Ut^2);
                d2Utp_dY2(1,index) = (Ut - ((Y - Yc)/(2*Ut)) )/(Ut^2);
                d2Utp_dZ2(1,index) = (Ut - ((Z - Zc)/(2*Ut)) )/(Ut^2);

                d2Utp_dXY(1,index) = -( (X - Xc)*(Y - Yc) ) / (Ut^3);
                d2Utp_dXZ(1,index) = -( (X - Xc)*(Z - Zc) ) / (Ut^3);
                d2Utp_dYZ(1,index) = -( (Y - Yc)*(Z - Zc) ) / (Ut^3);
            end
            
            alpha = d2Utp_dX2 * A;
            beta =  d2Utp_dY2 * A;
            gamma = d2Utp_dZ2 * A;
            
            xi =    d2Utp_dXY*A;
            phi =   d2Utp_dXZ*A;
            psi =   d2Utp_dYZ*A;
            
            E = (alpha.*alpha) + (beta.*beta) + (gamma.*gamma) + ...
                (2*xi.*xi) + (2*phi.*phi) + (2*psi.*psi);
             
            
            if(isnan(norm(E)))
                E = zeros(1,4);
                %warning('Bending energy is zero. This may happen for vertices lying on the seam of the planar surface');
           end
        end
        
        %%
        function [E, Phi_rr, Phi_tt, Phi_rt] = get_bending_energy(obj, A, D, control_handles, r, theta, N)
           theta = degtorad(theta);

           cT = cos(theta);
           sT = sin(theta);
           
           Phi_rr = zeros(1, 4, 'double');   
           Phi_tt = zeros(1, 4, 'double');
           Phi_rt = zeros(1, 4, 'double');
           
           Phi_rr_prio = zeros(1, N, 'double');   
           Phi_rt_prio = zeros(1, N, 'double');  
           Phi_tt_prio_1 = zeros(1, 4, 'double');  
           Phi_tt_prio_2 = zeros(1, N, 'double');  
           
           a1 = D(1,1); a5 = D(1,2); a9  = D(1,3);
           a2 = D(2,1); a6 = D(2,2); a10 = D(2,3);
           a3 = D(3,1); a7 = D(3,2); a11 = D(3,3);
           a4 = D(4,1); a8 = D(4,2); a12 = D(4,3);
           
           for j = 1:N
               xC = control_handles(j,1);
               yC = control_handles(j,2);
               zC = control_handles(j,3);
               
               U_t = obj.U( xC - sT, yC - r, zC - cT );
               
               
%                Phi_rr_prio(1,j) =   ((U_t).^2 - 2*( yC - r)) / (U_t).^3;
               Phi_rr_prio(1,j) =   ((U_t).^2 - ( r - yC).^2) / (U_t).^3;
               
               Phi_tt_prio_2(1,j) = ((U_t).^2*((zC*cT) + (xC*sT)) - (((zC*sT) - (xC*cT)).^2))/(U_t).^3;
               
%                Phi_rt_prio(1,j) =   (((zC*sT - xC*cT) - (r - yC))*(r - yC)) / (U_t).^(3/2);
               Phi_rt_prio(1,j) =  - (((zC*sT - xC*cT))*(r - yC)) / (U_t).^(3);
           end
           
           
           Phi_tt_prio_1(1,1) = -((a1*sT) + (a3*cT));
           Phi_tt_prio_1(1,2) = -((a5*sT) + (a7*cT));
           Phi_tt_prio_1(1,3) = -((a9*sT) + (a11*cT));
           
           Phi_rr = Phi_rr_prio * A;
           Phi_rt = Phi_rt_prio * A;
           %display(size(Phi_rt));
           Phi_tt = Phi_tt_prio_1 + (Phi_tt_prio_2 * A);
           
           E = Phi_rr*Phi_rr.' + (2*(Phi_rt*Phi_rt.')) + (Phi_tt*Phi_tt.');  
           
           if(isnan(E))
                E = 0;
               % warning('Bending energy is zero. This may happen for vertices lying on the seam of the planar surface');
           end
           
           %display('Dim of E');
           %display(size(E));
           %display((Phi_tt).^2);
           
        end
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Compute the 3D TPS parameters corresponding to target positions        % 
        %  for the control handles                                                %
        %                                                                         %
        % @param: control_handles - the input control handles for the TPS         %
        % @param: control_handle_target - the target position for the control     %
        %                                   handles of the TPS                    %
        % @param: N - number of 3D points in data                                 %
        %                                                                         % 
        % @return: [A, d] - computed parameters of the TPS                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [A, d, Linv, xi] = compute_tps(obj, control_handles, control_handle_target, N)
            K = zeros(N, N);
            P = ones (N, 4);
            O = zeros(4, 4);
            Y = zeros((N+4), 4);

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
            
%             fprintf('Rank: %d, rows %d\n',rank(K),N);
            
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
            
%             Uncomment this (VERY IMPORTANT!)
            E = Y(1:N,1).'*xi*Y(1:N,1) + Y(1:N,2).'*xi*Y(1:N,2) + Y(1:N,3).'*xi*Y(1:N,3);
            J = zeros(1,(3*N));
            
%             E = 0; %%% ----> DELETE this!
            
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
            
%             Uncomment this (VERY IMPORTANT!)
            J(1,:) = ( dYxT_dY * kron((xi*Y(1:N,1)), ones_block_thrice) + Y(1:N,1).'*xi*dYx_dY + ...
                        dYyT_dY * kron((xi*Y(1:N,2)), ones_block_thrice) + Y(1:N,2).'*xi*dYy_dY + ...
                        dYzT_dY * kron((xi*Y(1:N,3)), ones_block_thrice) + Y(1:N,3).'*xi*dYz_dY );
            
        end


        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Find the r, theta coordinates of a point lying on the right circular   %
        %  cylinder, assuming their principle axis aligns with the Y-axis and     %
        %  passes through the origin                                              % 
        %                                                                         %
        % @param: X, Y, Z - cartesian coordinates                                 %
        %                                                                         % 
        % @return: [r, theta] - computed values                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         function [r, theta] = cylindrical_coordinates_Yaligned_forward(obj, X, Y, Z)
% %             r = Y;
% %             dot_pr = (0*X) + (Z*1);
% %             norm_u = sqrt(1);
% %             norm_v = sqrt((X*X) + (Z*Z));
% %             CosTheta = dot_pr/(norm_u*norm_v);
% %             theta = acosd(CosTheta);    
% %             if(X < 0)
% %                 theta = 360.0 -theta;
% %             end
% %         end

        function [r, theta] = flattened_coordinates(obj, X, Y, Z)
            r = Y;
            dot_pr = (0*X) + (Z*1);
            norm_u = sqrt(1);
            norm_v = sqrt((X*X) + (Z*Z));
            CosTheta = dot_pr/(norm_u*norm_v);
            theta = acosd(CosTheta);    
            if(X < 0)
                theta = 360.0 -theta;
            end
%             fprintf('%d %d %d\n', X, Z, theta);
%             if(X > 0)
%                 theta = 360.0 -theta;
%             end
%             theta = -theta;
%             theta = deg2rad(theta);
        end


        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Find the normal to the parameterized surface                           %
        %                                                                         %
        % @param: A - A matrix corresponding to the TPS                           %
        % @param: D - D matrix corresponding to the TPS                           %
        % @param: control_handles - the input control handles for the TPS         %
        % @param: r - the r coordinate of the point for which the normal needs to % 
        %             be computed                                                 %
        % @param: theta - the theta coordinate of the point for which the         % 
        %              normal needs to be computed                                %
        % @param: N - number of 3D points in data                                 %
        %                                                                         % 
        % @return: [Nx, Ny, Nz] - the computed normals                            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Nx, Ny, Nz, lambda_r_prio, lambda_theta_prio, phi_r, phi_theta, norm, cT, sT] = get_normal(obj, A, D, control_handles, r, theta, N)
           theta = degtorad(theta);

           lambda_r_prio = zeros(1, N, 'double');    
           cT = cos(theta);
           sT = sin(theta);
           for j = 1:N
               lambda_r_prio(1,j) = (r - control_handles(j,2)) / ( obj.U( control_handles(j,1) -   sT,     ...
                                  control_handles(j,2) -   r,     ...
                                  control_handles(j,3) -   cT ));        
           end
            lambda_r = lambda_r_prio * A;
            lambda_1_r = lambda_r(1,1);
            lambda_2_r = lambda_r(1,2);
            lambda_3_r = lambda_r(1,3);
            lambda_4_r = lambda_r(1,4);

           lambda_theta_prio = zeros(1, N, 'double');    
           
           for j = 1:N
               lambda_theta_prio(1,j) = ((control_handles(j,3)*sT) - (control_handles(j,1)*cT)) / ( obj.U(control_handles(j,1) -   sT,     ...
                                  control_handles(j,2) -   r,     ...
                                  control_handles(j,3) -   cT ));
           end
           
            lambda_theta   = lambda_theta_prio * A;
            lambda_1_theta = lambda_theta(1,1);
            lambda_2_theta = lambda_theta(1,2);
            lambda_3_theta = lambda_theta(1,3); 
            lambda_4_theta = lambda_theta(1,4); 

            a1 = D(1,1); a5 = D(1,2); a9  = D(1,3);  a13 = D(1,4);
            a2 = D(2,1); a6 = D(2,2); a10 = D(2,3); a14 = D(2,4);
            a3 = D(3,1); a7 = D(3,2); a11 = D(3,3); a15 = D(3,4);
            a4 = D(4,1); a8 = D(4,2); a12 = D(4,3); a16 = D(4,4);
            
            theta_1_1 = (a1*cT) - (a3*sT);
            theta_1_2 = (a5*cT) - (a7*sT);
            theta_1_3 = (a9*cT) - (a11*sT);
            theta_1_4 = (a13*cT) - (a15*sT);
            
            phi_r = [(a2 + lambda_1_r) (a6 + lambda_2_r) (a10 + lambda_3_r) (a14 + lambda_4_r)];
            phi_theta = [(theta_1_1 + lambda_1_theta) (theta_1_2 + lambda_2_theta) (theta_1_3 + lambda_3_theta) (theta_1_4 + lambda_4_theta)];

            Nx = ((a6 + lambda_2_r)*((a9*cT) - (a11*sT) + lambda_3_theta)) - ...
                    ((a10 + lambda_3_r)*((a5*cT) - (a7*sT) + lambda_2_theta));
            Ny = -((a2 + lambda_1_r)*((a9*cT) - (a11*sT) + lambda_3_theta)) + ...
                    ((a10 + lambda_3_r)*((a1*cT) - (a3*sT) + lambda_1_theta));
            Nz = ((a2 + lambda_1_r)*((a5*cT) - (a7*sT) + lambda_2_theta)) - ...
                    ((a6 + lambda_2_r)*((a1*cT) - (a3*sT) + lambda_1_theta));    

            norm = sqrt((Nx*Nx) + (Ny*Ny) + (Nz*Nz));

            Nx = Nx/norm;
            Ny = Ny/norm;
            Nz = Nz/norm;   
        end

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  U value for the TPS                                                    %
        %  CAUTION: This is just a L2 norm function, and not the complete U value %
        %  The call to this function needs to handle generating the U value
        %  
        % @param: x,y, z - cartesian coordinates of the point                     %
        %                                                                         %
        % @return: Ur - computed U-value for this point                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Ur = U(obj, x,y,z)
            Ur = sqrt((x*x) + (y*y) + (z*z));
            if(abs(Ur) < 0.0000001)
                Ur = 0.0001;
            end
        end

        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Find correspondence between arbitrary 3D point and cylindrical template%
        %  using direct projection to central/principal axis of cylinder,
        %  assuming the cylinder is aligned with the Y-axis.
        %                                                                         %
        % @param: target - [Nx3] matrix containing a list of 3D points in         %   
        % cartesian coordinates                                                   %
        %                                                                         %
        % @return: [Nx3] matrix containing cartesian coordinates on unit, circular% 
        %                cylinder                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function target_natural = find_correspondence_natural(obj, target, N, radius)
            if ~exist('radius','var')
                radius = 1.0;
            end
            target_natural = zeros(N, 3, 'double');  
            for i = 1:N
                X = target(i,1); Z = target(i,3);
                eucl_dist = sqrt((X*X) + (Z*Z));
                X = (radius*X) / eucl_dist;
                Z = (radius*Z) / eucl_dist;
                target_natural(i,1) = X;
                target_natural(i,2) = target(i,2);
                target_natural(i,3) = Z;
            end
        end
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Find correspondence between arbitrary 3D point and a rectangular lattice grid%
        %  using a fixed number of U, V, W resolution
        %  such that the entire data volume is encompassed in the lattice
        %  grid
        % @param: target - [Nx3] matrix containing a list of 3D points in         %   
        % cartesian coordinates                                                   %
        % @param: N: number of 3D points in 'target'
        % @param: [n_pts_X, n_pts_Y, n_pts_Z]: number of partitions to the
        % lattice grid in X, Y, Z direction (corresponding to the U, V, W
        % terms in Blender)
        % @param: padding: padding to the bounding box
        %                                                                         %
        % @return: [Nx3] matrix containing cartesian coordinates on unit, circular% 
        %                cylinder                                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function target_natural = find_correspondence_lattice_grid(obj, target, N, n_pts_X, n_pts_Y, n_pts_Z, padding)
            if ~exist('n_pts_X', 'var')
                n_pts_X = 5;
            end
            if ~exist('n_pts_Y', 'var')
                n_pts_Y = 10;
            end
            if ~exist('n_pts_Z', 'var')
                n_pts_Z = 5;
            end     
            if ~exist('padding', 'var')
                padding = 0.1;
            end     
            n_pts_X = n_pts_X - 1;n_pts_Y = n_pts_Y - 1;n_pts_Z = n_pts_Z - 1;
            target_natural = zeros((n_pts_X * n_pts_Y * n_pts_Z), 3, 'double');  
            x_min = realmax;
            y_min = realmax;
            z_min = realmax;
            
            x_max = -realmax;
            y_max = -realmax;
            z_max = -realmax;
            
            for i = 1:N
                X = target(i,1);
                Y = target(i,2);
                Z = target(i,3);
                if(X > x_max)
                    x_max = X;
                end
                if(X < x_min)
                    x_min = X;
                end
                
                if(Y > y_max)
                    y_max = Y;
                end
                if(Y < y_min)
                    y_min = Y;
                end
                
                if(Z > z_max)
                    z_max = Z;
                end
                if(Z < z_min)
                    z_min = Z;
                end
            end
            
            x_max = x_max + padding;
            y_max = y_max + padding;
            z_max = z_max + padding;
            x_min = x_min - padding;
            y_min = y_min - padding;
            z_min = z_min - padding;            
            
            X = x_min;
            Y = y_min;
            Z = z_min;
            
            x_step = (x_max - x_min)/n_pts_X;
            y_step = (y_max - y_min)/n_pts_Y;
            z_step = (z_max - z_min)/n_pts_Z;
            index = 1;
            
            for index1 = 1:(n_pts_X+1)
                for index2 = 1:(n_pts_Y+1)
                    for index3 = 1:(n_pts_Z+1)
                        target_natural(index,1) = X;
                        target_natural(index,2) = Y;
                        target_natural(index,3) = Z;
                        Z = Z + z_step;
                        index = index + 1;
                    end
                    Y = Y + y_step;
                    Z = z_min;
                end
                X = X + x_step;
                Y = y_min;
                Z = z_min;
            end
            
        end
        
        
        %%
%         Factorises the input (N) into three integers, the product of which equals to the input 'N' (if possible)
        function [num_X, num_Y, num_Z] = factorise_pointcloud_size(obj, N)
            
            F = factor(N);
            dim_f = size(F);
            dim_f = dim_f(2);
            if(dim_f < 3)
                fprintf('No. of vertices in input mesh cannot be factored into 3 dimensions with ineteger multiplicand, this is going to produce incorrect results. Please change the no. of vertices in input mesh and retry');
            end
            
            if(dim_f <= 4)
                num_X = F(1);
                num_Y = 1.0;
                num_Z = F(2);
                for nI = 3:dim_f
                    num_Y = num_Y * F(nI);
                end
            else
                num_X = F(1)*F(2);
                num_Y = 1.0;
                num_Z = F(3)*F(4);
                for nI = 5:dim_f
                    num_Y = num_Y * F(nI);
                end
            end
        end
        


        %% EXPERIMENTS with a cylinder in 3D
        % Call from runme_exp1.m -> uncomment the last two lines (while
        % commenting out the rest)
        function example(obj)
            formatSpec = '%f %f %f';
            sizeData = [3 Inf];
            fileID_ip = fopen('data/additional/HighDensity.xyz','r'); 
            control_handles = fscanf(fileID_ip,formatSpec,sizeData);
            fclose(fileID_ip);
            fileID_target = fopen('data/SDP_Quadrics/test_piece_rot.xyz','r'); 
            control_handle_target = fscanf(fileID_target,formatSpec,sizeData);
            fclose(fileID_target);
            
            control_handles = control_handles';
            control_handle_target = control_handle_target';                        

            dim = size(control_handles);
            N = dim(1);

            dim_target = size(control_handle_target);
            N_target = dim_target(1);
            
            
%             [n_x, n_y, n_z] = obj.factorise_pointcloud_size(N_target);
            
            base = obj.find_correspondence_natural(control_handle_target, N_target);

            [A, d, Linv, xi] = obj.compute_tps(base, control_handle_target, N_target);             

            
%             [n_x, n_y, n_z] = obj.factorise_pointcloud_size(N);
            
            base2 = obj.find_correspondence_natural(control_handles, N);
                                    
            for index = 1:30                   
            
             %% OBSERVATIONS:
                %%% -> Moving control handle (control_handle_target) works for sure
                %%% moving the base (base2) works as well
                %%% but moving the other target (base) does NOT work at all
                %%% also: there is difference between moving control handle
                %%% vs base, i.e., control handle -> more obvious; base ->
                %%% subtle.
                %%% -> 'interpolate_tps' could be used for interpolating
                %%% high-resolution template, check method comments on top.
                %%% -> Preferable to control deformation using control handle
                %%% (control_handle_target) for proper normal estimation,
                %%% i.e., moving the control handle is better for normals
                %%% -> So the two pointcloud arguments of interpolate_tps are best kept
                %%% the same (later shown, not necesssarily)
             %%
%                 if(index<11)
%                     axis = 1;
%                     for s_index = 1:10
%                         control_handle_target(s_index,axis) = control_handle_target(s_index,axis) + 0.01;
%                     end
%                 elseif((index > 10) && (index < 21))
%                     axis = 2;
%                     for s_index = 1:10
%                         control_handle_target(s_index,axis) = control_handle_target(s_index,axis) + 0.01;
%                     end
%                 elseif((index > 20))
%                     axis = 3;
%                     for s_index = 1:10
%                         control_handle_target(s_index,axis) = control_handle_target(s_index,axis) + 0.01;
%                     end
%                 end                                        

%                 if(index<11)
%                     axis = 1;
%                     for s_index = 1:10
%                         base2(s_index,axis) = base2(s_index,axis) + 0.01;
%                     end
%                 elseif((index > 10) && (index < 21))
%                     axis = 2;
%                     for s_index = 1:10
%                         base2(s_index,axis) = base2(s_index,axis) + 0.01;
%                     end
%                 elseif((index > 20))
%                     axis = 3;
%                     for s_index = 1:10
%                         base2(s_index,axis) = base2(s_index,axis) + 0.01;
%                     end
%                 end       

%                      base2(15,3) = base2(15,3) + 0.02;
%                      control_handle_target(15,3) = control_handle_target(15,3) + 0.001;

                [A, d, E, J] = obj.update_tps(control_handle_target, Linv, xi, N_target);
                
                base2(20:2180,:) = [];
                N = obj.length_matrix(base2);
                    
                [Px, Py, Pz, Nx, Ny, Nz, E, jacobian_components, dQ_dxi, C] = obj.interpolate_tps_global(A, d, base2, ...
                                                                        base, N, N_target, xi, Linv, 1);
                                                                    
                data_post = horzcat(Px, Py, Pz, Nx, Ny, Nz);
                
                 [data_post] = obj.color_data_with_curvature(data_post, C);
                                                                    
                final_path = strcat('/home/agniva/Documents/figures/Adrien_W8.21/TPS_Interpolation/test_',num2str(index),'.xyz');

                fileID_OP = fopen(final_path,'w');                                                                
                for i = 1:N
                   fprintf(fileID_OP,'%5d %5d %5d %5d %5d %5d\n', Px(i,1), Py(i,1), Pz(i,1), ...
                                                                 Nx(i,1), Ny(i,1), Nz(i,1) );
                end
                fclose(fileID_OP);
                
                base_path = strcat('/home/agniva/Documents/figures/Adrien_W8.21/TPS_Interpolation/base_',num2str(index),'.xyz');

                fileID_OP = fopen(base_path,'w');                                                                
                for i = 1:N_target
                   fprintf(fileID_OP,'%5d %5d %5d\n',control_handle_target(i,1), control_handle_target(i,2), control_handle_target(i,3));
                end
                fclose(fileID_OP);
                
            end

        end
        
        function [rotated_data] = color_data_with_curvature(obj, data, curvature)

            should_paint_curvature = 1;
            if ~exist('curvature','var')
                should_paint_curvature = 0;
            end
%             for frame = 1:N_frames
                size = obj.length_matrix(data);
                
                if(should_paint_curvature == 1)
                    C = abs(curvature);
                    C(C<0.19)  = 0; % thresholding out small curvatures
                    m = min(C(:));
                    C = (C(:)-m)/(max(C(:))-m)*10000;
                end

                for index = 1: size
                    P = zeros(3,1);
                    P(1,1) = data(index,1);
                    P(2,1) = data(index,2);
                    P(3,1) = data(index,3);
                    rotated_data_temp(index,1) = P(1,1);
                    rotated_data_temp(index,2) = P(2,1);
                    rotated_data_temp(index,3) = P(3,1);

                    if(should_paint_curvature == 1)
                        N = zeros(3,1);
                        N(1,1) = data(index,4);
                        N(2,1) = data(index,5);
                        N(3,1) = data(index,6);
                        rotated_data_temp(index,4) = N(1,1);
                        rotated_data_temp(index,5) = N(2,1);
                        rotated_data_temp(index,6) = N(3,1);           
                        if(C(index) > 255)
                            val = 255;
                        else
                            val = C(index);
                        end
                        rotated_data_temp(index,7) = val;
                        rotated_data_temp(index,8) = 0;%255.0 - C(index);
                        rotated_data_temp(index,9) = 0;%C(index);
                        rotated_data_temp(index,10) = 1.0;
                    end
                end
                rotated_data = rotated_data_temp;
%             end

        end
        
        function [N] = length_matrix(obj, X)
            N = size(X);
            N = N(1);
        end

    end
end

