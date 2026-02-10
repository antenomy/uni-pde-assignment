clear; 
close all; 
clc;

epsilon = 0.1;
N = 250;
max_iterations = 2000000;
tol = 1e-10;

deltas = [1e-4, 1e-5];

figure; hold on; grid on;
xlabel('iterations'); ylabel('log_{10}(step norm)');
title(sprintf('log_{10} step norm (\\epsilon=%.2f, N=%d)', epsilon, N));

for k = 1:numel(deltas)
    delta = deltas(k);
    [X, u, iters, step_history, deltaBound] = solve_tend_gradient(epsilon, N, delta, max_iterations, tol);
    
    sh = step_history(:);
    n = numel(sh);
    finiteMask = isfinite(sh) & (sh > 0);
    idx = find(finiteMask);
    if ~isempty(idx)
        plot(idx, log10(sh(idx)), 'DisplayName', sprintf('\\delta = %.0e', delta));
    end
    
end

legend('Location','best');
hold off;

saveas(gcf, 'p7_log10stepnorm.png');

%[X, u, iters, step_history, deltaBound] = solve_tend_gradient(epsilon, N, delta, max_iterations, tol);
%figure;
%hold on; 
%grid on;
%plot(X, u); 
%xlabel('r');
%ylabel('U(r)');
%hold off;

%figure;
%hold on; 
%grid on;
%semilogy(step_history); 
%xlabel('iterations');
%ylabel('step norm');
%hold off;

%disp(iters);
%disp(step_history(end));
%disp(deltaBound);



function [X, u, iters, step_history, deltaBound] = solve_tend_gradient(epsilon, N, delta, max_iterations, tol)
    [nodes_ref, weights] = gauss_ref();

    num_edges = N + 1;
    X = get_mesh_nodes(epsilon, num_edges);
    h = (1 - epsilon) / N;

    dphi_dx   = dphi_dx_from_h(h, nodes_ref);
    nodes_act = gauss_act(X, nodes_ref, N);

    stiff_local  = compute_stiff_local(h, dphi_dx, nodes_act, weights, N);
    stiff_global = compute_stiff_global(stiff_local, N);
    A = sparse(stiff_global);

    u = zeros(num_edges, 1);
    u(1) = 1;
    u(end) = 0;

    index = 2:N;
    b = -A(index, 1)*u(1) - A(index, end)*u(end);

    % Gradient method
    step_history = zeros(max_iterations, 1);
    A_int = A(index, index); % A interior, keep boundary fixed
    
    deltaBound = 1 / eigs(A_int, 1, 'largestreal');

    for n = 1:max_iterations
        gradient = 2 * (A_int * u(index) - b);
        u_next = u;
        u_next(index) = u(index) - delta * gradient;
        step_history(n) = norm(u_next(index) - u(index), 2);
        u = u_next;

        if step_history(n) < tol
            iters = n;
            step_history = step_history(1:n); % trim list
            return;
        end
    end
    iters = max_iterations;
end

function X = get_mesh_nodes(epsilon, num_edges)
    X = linspace(epsilon, 1, num_edges)';
end

function [gauss_nodes_ref, gauss_weights] = gauss_ref()
    % Just one point for this one
    gauss_nodes_ref = 0.5;
    gauss_weights = 1.0;
end

function gauss_nodes_act = gauss_act(X, gauss_nodes_ref, N)
    num_nodes = numel(gauss_nodes_ref);
    gauss_nodes_act = zeros(N, num_nodes);

    for i = 1:N 
        xL = X(i);
        xR = X(i+1);
        for w = 1:num_nodes
            gauss_nodes_act(i, w) = xL + (xR-xL)*gauss_nodes_ref(w);
        end
    end
end

function a_val = a_def(x)
    a_val = x; 
end

function dphi_dx = dphi_dx_from_h(h, gauss_nodes_ref)
    num_nodes = numel(gauss_nodes_ref);
    dphi_dx = repmat([-1; 1] / h, 1, num_nodes);  % 2 x num_nodes
end

function stiff_local = compute_stiff_local(h, dphi_dx, gauss_nodes_act, weights, N)
    num_nodes = numel(weights);
    stiff_local = zeros(N, 4);

    for i = 1:N
        A11 = 0; 
        A12 = 0;
        A21 = 0;
        A22 = 0;
        for w = 1:num_nodes
            node = gauss_nodes_act(i, w);

            coeff = weights(w) * a_def(node) * h;
            d1 = dphi_dx(1, w);
            d2 = dphi_dx(2, w);

            A11 = A11 + coeff * d1 * d1;
            A12 = A12 + coeff * d1 * d2;
            A21 = A21 + coeff * d2 * d1;
            A22 = A22 + coeff * d2 * d2;
        end
        
        stiff_local(i, :) = [A11 A12 A21 A22];
    end
end

function stiff_global = compute_stiff_global(stiff_local, N)
    stiff_global = zeros(N + 1);

    for i = 1:N
        j = i + 1;

        stiff_global(i, i) = stiff_global(i, i) + stiff_local(i, 1); %A11
        stiff_global(i, j) = stiff_global(i, j) + stiff_local(i, 2); %A12
        stiff_global(j, i) = stiff_global(j, i) + stiff_local(i, 3); %A21
        stiff_global(j, j) = stiff_global(j, j) + stiff_local(i, 4); %A22
    end

    
end