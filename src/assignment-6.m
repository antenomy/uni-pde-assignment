epsilons = [0.5, 0.2, 0.1, 0.05, 0.02];
N=250;


figure;
hold on;
grid on;

for m = 1:numel(epsilons)
    [X, u] = solve_tent_fem(epsilons(m), N);
    plot(X, u, '-');
    %solve_tent_fem(epsilon, N)

    %epsilon = epsilons(m);
    %h = (1-epsilon)/N;
    
    %dof_pos
end

legend(string(epsilons));
xlabel('r');
ylabel('U(r)');



function [X, u] = solve_tent_fem(epsilon, N)
    % Gauss nodes and weights
    [nodes_ref, weights] = gauss_ref();

    num_edges = N + 1;
    % mesh
    X = get_mesh_nodes(epsilon, num_edges);
    h = (1-epsilon) / N;
    
    % We use hat functions \phi for our basis, d
    % For each "element" i.e. interval [x_{i-1}, x_i], we gotta consider
    % the left hat function sloping down and right one rising up. We encode
    % this here:
    dphi_dx = dphi_dx_from_h(h, nodes_ref);

    nodes_act = gauss_act(X, nodes_ref, N);
    
    % stifflocal is N x 4
    stiff_local = compute_stiff_local(h, dphi_dx, nodes_act, weights, N);
    
    stiff_global = compute_stiff_global(stiff_local, N);
    A = sparse(stiff_global);

    % No load vector (RHS) for this problem
    u = zeros(num_edges, 1);
    u(1) = 1;
    u(end) = 0;

    index = 2:N;
    b = -A(index, 1)*u(1) - A(index, end)*u(end);
    u(index) = A(index, index) \ b;
    
    %b = zeros(num_edges, 1);
    %for i = 2:N
    %    b(i) = -A(i, 1)*u(1) - A(i, end)*u(end);
    %    u(i) = A \ b;
    %end


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