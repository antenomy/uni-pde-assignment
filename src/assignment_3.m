N_values = [10, 20, 40, 80, 160, 320];

g = 1;
% comparison value
a_func = @(x) 1 + x; 
f_func = @(x) 0;
u_exact = @(x) log(1 + x);

p = zeros(2, 1);

% a_func = @(x) exp(x); 
% f_func = @(x) exp(x);
% u_exact = @(x) -(1+exp(1)).*exp(-x) +1 +exp(1) - x;

figure(1);
for k = 1:2
    errors = zeros(size(N_values));
    for i = 1:length(N_values)
        N = N_values(i);
        
        [X, U] = fem_1d_heat(a_func, f_func, g, N, k);
        
        % calculate exact U at these node points
        U_true = u_exact(X);
        
        % calculate L2 Error Norm
        diff = U - U_true;
        dx = 1/N;
        
        errors(i) = sqrt(sum(diff.^2) * dx);
    end

    loglog(N_values, errors, '-s');
    p_temp = polyfit(log(1./N_values), log(errors), 1)
    p(k) = p_temp(1);
    hold on;
end

grid on;
legend(["1-point, slope " + p(1), "2-point, slope " + p(2)]);
xlabel('Number of Elements (N)');
ylabel('L2 Error Norm');
title('Convergence Rate Analysis');
hold off;


function [X, U, K, F] = fem_1d_heat(a_func, f_func, g, N, point_count)
%   FEM_1D_HEAT Solves a 1D steady-state heat equation using the Finite Element Method.
%
%   INPUTS:
%       a_func      - Function handle for thermal conductivity a(x).
%       f_func      - Function handle for heat source term f(x).
%       g           - Scalar value for the Neumann boundary condition at x=1.
%       N           - Integer number of elements (intervals) to use.
%       point_count - Integer (1 or 2) specifying the Gauss quadrature rule 
%                     to use for numerical integration.
%
%   OUTPUTS:
%       X           - (N+1)x1 vector of nodal coordinates.
%       U           - (N+1)x1 vector of the solution u(x) at each node.
%       K           - NxN Global Stiffness Matrix (after removing fixed DOF).
%       F           - Nx1 Global Load Vector (after removing fixed DOF).
%

    X = linspace(0, 1, N+1)';
    dx = 1/N;
    dphi_dx = [-1/dx; 1/dx];
    
    % global stiffness & load
    K = zeros(N, N);
    F = zeros(N, 1);

    % define Gauss points and weights
    if point_count == 2
        xi_list = [-1/sqrt(3), 1/sqrt(3)];
        w_list  = [1, 1]; 
    else
        xi_list = [0];
        w_list  = [2];
    end

    for i = 1:N
        x_i = X(i); 

        % local stiffness & load 
        k_local = zeros(2,2);
        f_local = zeros(2,1);

        % loop over Gauss points
        for q = 1:length(xi_list)
            xi = xi_list(q);
            w  = w_list(q);

            x_q = (x_i + dx/2) + (dx/2) * xi;

            phi = [ (1 - xi)/2; 
                    (1 + xi)/2 ];

            a_val = a_func(x_q);
            f_val = f_func(x_q);

            
            k_local = k_local + w * a_val * (dphi_dx * dphi_dx');
            f_local = f_local + w * f_val * phi;
        end

        k_local = k_local * (dx/2);
        f_local = f_local * (dx/2);
        nodes = [i-1, i];
        
        % add all integrated pieces into appropriate place in global stiffness/load matrices
        for p = 1:2
            i_p = nodes(p);
            if i_p ~= 0
                F(i_p) = F(i_p) + f_local(p);
                for q = 1:2
                    i_q = nodes(q);
                    if i_q ~= 0
                        K(i_p, i_q) = K(i_p, i_q) + k_local(p, q);
                    end
                end
            end
        end
    end
    
    F(end) = F(end) + g; %Neumann BC 

    U_inner = K \ F;
    
    U = zeros(N+1, 1);
    U(1) = 0; % Dirichlet BC
    U(2:end) = U_inner;
end