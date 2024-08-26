% Parameters
P = 50;
mu = 1;
sigma = 4;
N = 100; % Number of elements in the finite element grid
h = 1/N; % Element size
M = 1000; % Number of Monte Carlo samples

% Finite element grid
x = linspace(0, 1, N+1)'; % N+1 nodes in the interval [0, 1]

% Assembly of the stiffness matrix and load vector
A = zeros(N+1, N+1); % Stiffness matrix
b = h * ones(N+1, 1); % Load vector for f = 1

% Basis function derivatives (constant in each element)
phi_prime = [-1/h; 1/h];

% Assembling the stiffness matrix A
for i = 1:N
    A(i:i+1, i:i+1) = A(i:i+1, i:i+1) + h * (phi_prime * phi_prime');
end

% Apply boundary conditions
A(1, :) = 0; A(1, 1) = 1; % u(0) = 0
A(end, :) = 0; A(end, end) = 1; % u(1) = 0
b(1) = 0; b(end) = 0;

% Monte Carlo Simulation
Q_values = zeros(M, 1);

for m = 1:M
    % Generate random variables xi_p(w)
    xi = -1 + 2 * rand(P, 1); % Uniform distribution U(-1, 1)
    
    % Compute a(x, omega) on the grid
    a = mu + sum(sigma * (cos((1:P)' * pi * x') ./ ((1:P)'.^2 * pi^2)) .* xi, 1)';
    
    % Modify the stiffness matrix with a(x, omega)
    A_mod = A;
    for i = 1:N
        A_mod(i:i+1, i:i+1) = A_mod(i:i+1, i:i+1) .* a(i);
    end
    
    % Solve the linear system for u_h(x, omega)
    u_h = A_mod \ b;
    
    % Compute the integral mean Q(u_h)
    Q_values(m) = mean(u_h);
end

% Estimate the expected value E[Q(u)]
E_Q_u = mean(Q_values);

% Display the result
fprintf('Estimated E[Q(u)] with M = %d samples: %f\n', M, E_Q_u);

% Convergence study (optional)
% Increase M and/or decrease h to study the convergence of E[Q(u)]
