% Parameters
P = 50;
mu = 1;
sigma = 4;
N = 100; % Number of elements in the finite element grid
h = 1/N; % Element size

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

% Generate 20 random samples and compute solutions
num_samples = 20;
solutions = zeros(N+1, num_samples);

for m = 1:num_samples
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
    
    % Store the solution
    solutions(:, m) = u_h;
end

% Plot the solutions
figure;
hold on;
for m = 1:num_samples
    plot(x, solutions(:, m), 'DisplayName', ['\omega_' num2str(m)]);
end
xlabel('x');
ylabel('u_h(x, \omega_m)');
title('Approximate Solutions for 20 Random Samples');
legend;
hold off;
