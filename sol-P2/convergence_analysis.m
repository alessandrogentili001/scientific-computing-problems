% Parameters
P = 50;
mu = 1;
sigma = 4;

% Convergence study parameters
M_values = [100, 500, 1000, 5000, 10000]; % Monte Carlo samples
N_values = [50, 100, 200, 400, 800]; % Number of elements for mesh refinement

% Storage for results
E_Q_u_M = zeros(length(M_values), 1);
E_Q_u_N = zeros(length(N_values), 1);

% Convergence with respect to M (fixed N)
N = 100; % Fixed mesh size
h = 1/N; % Element size

x = linspace(0, 1, N+1)'; % N+1 nodes in the interval [0, 1]
A = zeros(N+1, N+1); % Stiffness matrix
b = h * ones(N+1, 1); % Load vector for f = 1
phi_prime = [-1/h; 1/h];

for i = 1:N
    A(i:i+1, i:i+1) = A(i:i+1, i:i+1) + h * (phi_prime * phi_prime');
end

A(1, :) = 0; A(1, 1) = 1; % u(0) = 0
A(end, :) = 0; A(end, end) = 1; % u(1) = 0
b(1) = 0; b(end) = 0;

for k = 1:length(M_values)
    M = M_values(k);
    Q_values = zeros(M, 1);

    for m = 1:M
        xi = -1 + 2 * rand(P, 1); % Uniform distribution U(-1, 1)
        a = mu + sum(sigma * (cos((1:P)' * pi * x') ./ ((1:P)'.^2 * pi^2)) .* xi, 1)';
        A_mod = A;
        for i = 1:N
            A_mod(i:i+1, i:i+1) = A_mod(i:i+1, i:i+1) .* a(i);
        end
        u_h = A_mod \ b;
        Q_values(m) = mean(u_h);
    end

    E_Q_u_M(k) = mean(Q_values);
end

% Convergence with respect to N (fixed M)
M = 1000; % Fixed number of Monte Carlo samples

for k = 1:length(N_values)
    N = N_values(k);
    h = 1/N; % Element size
    x = linspace(0, 1, N+1)'; % N+1 nodes in the interval [0, 1]
    A = zeros(N+1, N+1); % Stiffness matrix
    b = h * ones(N+1, 1); % Load vector for f = 1
    phi_prime = [-1/h; 1/h];
    
    for i = 1:N
        A(i:i+1, i:i+1) = A(i:i+1, i:i+1) + h * (phi_prime * phi_prime');
    end
    
    A(1, :) = 0; A(1, 1) = 1; % u(0) = 0
    A(end, :) = 0; A(end, end) = 1; % u(1) = 0
    b(1) = 0; b(end) = 0;
    
    Q_values = zeros(M, 1);

    for m = 1:M
        xi = -1 + 2 * rand(P, 1); % Uniform distribution U(-1, 1)
        a = mu + sum(sigma * (cos((1:P)' * pi * x') ./ ((1:P)'.^2 * pi^2)) .* xi, 1)';
        A_mod = A;
        for i = 1:N
            A_mod(i:i+1, i:i+1) = A_mod(i:i+1, i:i+1) .* a(i);
        end
        u_h = A_mod \ b;
        Q_values(m) = mean(u_h);
    end

    E_Q_u_N(k) = mean(Q_values);
end

% Plotting results

% Convergence with respect to M
figure;
plot(M_values, E_Q_u_M, '-o');
xlabel('Number of Monte Carlo samples (M)');
ylabel('Estimated E[Q(u)]');
title('Convergence of E[Q(u)] with respect to M');
grid on;

% Convergence with respect to N (mesh size)
figure;
plot(1./N_values, E_Q_u_N, '-o');
xlabel('Mesh size (1/h)');
ylabel('Estimated E[Q(u)]');
title('Convergence of E[Q(u)] with respect to Mesh Size (h)');
grid on;
