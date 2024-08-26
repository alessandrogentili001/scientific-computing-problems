function helmholtz_fem_solver
    % parameters 
    k = 2*pi;  % Wave number
    N_values = [10, 20, 40, 80, 160, 320];  % Number of elements to test

    % initialize errors arrays 
    errors_L2 = zeros(size(N_values));
    errors_H1 = zeros(size(N_values));

    figure('Position', [100, 100, 1200, 400]);

    for i = 1:length(N_values)
        N = N_values(i);
        [x, u_h, u_exact, err_L2, err_H1] = solve_helmholtz(k, N);

        errors_L2(i) = err_L2;
        errors_H1(i) = err_H1;

        % Plot solution comparison
        subplot(1, 2, 1);
        plot(x, real(u_h), 'b-', x, real(u_exact), 'r--', ...
             x, imag(u_h), 'g-', x, imag(u_exact), 'm--');
        hold on;
    end

    % Finalize solution comparison plot
    subplot(1, 2, 1);
    legend('Re(u_h)', 'Re(u_{exact})', 'Im(u_h)', 'Im(u_{exact})');
    xlabel('x');
    ylabel('u');
    title('FEM Solution vs Exact Solution');
    hold off;

    % Plot errors
    subplot(1, 2, 2);
    loglog(N_values, errors_L2, 'bo-', N_values, errors_H1, 'rs-');
    legend('L2 Error', 'H1 Error');
    xlabel('Number of Elements');
    ylabel('Error');
    title('Convergence of FEM Solution');
    grid on;

    % Display convergence rates
    L2_rate = polyfit(log(N_values), log(errors_L2), 1);
    H1_rate = polyfit(log(N_values), log(errors_H1), 1);
    fprintf('L2 convergence rate: %f\n', -L2_rate(1));
    fprintf('H1 convergence rate: %f\n', -H1_rate(1));
end