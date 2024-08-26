function [x, u_h, u_exact, err_L2, err_H1] = solve_helmholtz(k, N)
    h = 1/N;  % Element size
    x = linspace(0, 1, N+1)';

    % Assemble system matrices
    K = sparse(N+1, N+1);
    M = sparse(N+1, N+1);
    for e = 1:N
        nodes = [e, e+1];
        xe = x(nodes);
        Ke = local_stiffness(xe);
        Me = local_mass(xe);
        
        K(nodes, nodes) = K(nodes, nodes) + Ke;
        M(nodes, nodes) = M(nodes, nodes) + Me;
    end

    % Apply boundary conditions
    K(1,:) = 0;
    K(1,1) = 1;
    M(1,:) = 0;
    
    % Correct Robin boundary condition
    K(end,end) = K(end,end) + 1i*k;
    
    % Assemble system
    A = K - k^2*M;
    
    % Set up right-hand side
    b = zeros(N+1, 1);
    b(end) = k*exp(-1i*k);  % Contribution from Robin condition

    % Solve system
    u_h = A\b;

    % Compute exact solution
    u_exact = sin(k*x);

    % Compute errors
    err_L2 = sqrt(h * sum(abs(u_h - u_exact).^2));
    err_H1 = sqrt(err_L2^2 + h * sum(abs(diff(u_h)/h - k*cos(k*x(1:end-1))).^2));
end