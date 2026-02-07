

% -(a(x)u')' = f(x) on (0,1)
% BC: u(0)=0, a(1)u'(1)=g, g constant (I think)

N = 40;
g = 1;

% Specific a and f used in part 3
a_func = @(x) 1 + x; 
f_func = @(x) 0;

%a_func = @(x) exp(x);
%f_func = @(x) exp(x);

[X, U] = fem_1d_heat(a_func, f_func, g, N);
figure;
plot(X, U, '-o');

% Currently only does Gauss 1-point rule atm
function [X, U, K, F] = fem_heat_1d(a_func, f_func, g, N)
    X = linspace(0, 1, N+1)';
    dx = 1/N;
    
    %global stiffness & load
    K = zeros(N, N);
    F = zeros(N, 1);

    for i = 1:N
        x_i = X(i); 
        %xi1 = X(i+1);

        % local stiffness & load 
        k_local = zeros(2,2);
        f_local = zeros(2,1);

        dphi_dx = [-1/dx; 1/dx];

        % TODO: add 2-point Gauss rule
        int_a = 0;
        x_k = x_i + dx/2;
        phi = [1/2; 1/2];
        int_a = int_a + 2 * a_func(x_k);
        f_local = f_local + 2 * f_func(x_k) * phi;

        int_a = int_a * (dx/2);
        f_local = f_local * (dx/2);

        k_local = (dphi_dx * dphi_dx') * int_a;  
        nodes = [i-1, i];
        
        % Add all integrated pieces into appropriate place in 
        % global stiffness/load matrices
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

% N = 20;
% x = sort(rand(N-1, 1));
% x = [0;x;1];
% Adiag = 1./(x(3:end) - x(2:end-1)) + 1./(x(2:end-1) - x(1:end-2));
% Aoffd = 1./(x(3:end-1) - x(2:end-2));
% A = diag(Adiag) - diag(Aoffd, 1) - diag(Aoffd, -1);
% plot(x, 'o');
% 
% % f = 1
% b = (x(3:end) - x(1:end-2)) / 2;
% u = zeros(N+1, 1);
% u(2:end-1) = A \ b;
% plot(x, u , '-o');
