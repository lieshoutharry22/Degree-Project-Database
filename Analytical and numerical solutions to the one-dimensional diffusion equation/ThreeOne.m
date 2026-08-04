function ThreeOne
    % Define parameters
    Nx = 129; % Number of spatial grid points
    L = 2000; % Length of the mineshaft in meters
    dx = L / (Nx - 1); % Spatial step size
    D = 5e-1; % Diffusivity of CO in m^2/s
    dt = 60; % Time step in seconds
    t_final = 26.5 * 3600; % Final time in seconds (26.5 hours converted to seconds)
    
    % Generate the spatial grid
    x = linspace(0, L, Nx);
    
    % Calculate numerical solution without source term for comparison
    c_numerical = numericalSolutionFTCS(Nx, dt, dx, D, t_final);
    
    % Calculate numerical solution with source term
    c_numerical_ventilated = ventilated_numericalSolutionFTCS(Nx, dt, dx, D, t_final, x);
    
    % Plotting
    figure;
    plot(x, c_numerical, 'LineWidth', 2, 'DisplayName', 'Unventilated Mineshaft');
    hold on;
    plot(x, c_numerical_ventilated, 'LineWidth', 2, 'DisplayName', 'Ventilated Mineshaft');
    hold off;
    xlabel('Position along the mineshaft (m)');
    ylabel('CO Mass Fraction');
    title('CO Concentration Profile in the Mineshaft at T = 26.5 hours');
    legend('show');
    grid on;
end

function c_numerical = ventilated_numericalSolutionFTCS(Nx, dt, dx, D, t_final, x)
    % Initialize concentration profile
    c_numerical = zeros(Nx, 1); % Solution vector
    c_numerical(1) = 1; % Initial CO mass fraction fixed at 1.0 at x = 0
    
    % Calculate coefficient for the FTCS scheme
    sigma = D * dt / dx^2;
    
    % Discretize the source term f(x,t) = -0.001e^(-0.01 x^(1.5))
    f = -0.001 * exp(-0.01 * x.^1.5);
    
    % Time-stepping loop
    for t = dt:dt:t_final
        % Initialize temporary vector for the new time step
        c_new = c_numerical;
        
        % Loop over space
        for i = 2:Nx - 1
            % Apply FTCS scheme with source term
            c_new(i) = c_numerical(i) + sigma * (c_numerical(i - 1) - 2 * c_numerical(i) + c_numerical(i + 1)) + dt * f(i);
        end
        
        % Apply boundary conditions
        c_new(1) = 1; % CO mass fraction fixed at 1.0 at x = 0
        c_new(end) = 0; % CO mass fraction fixed at 0 at x = L
        
        % Update concentration profile
        c_numerical = c_new;
    end
end

function c_numerical = numericalSolutionFTCS(Nx, dt, dx, D, t_final)
    % Initialize concentration profile
    c_numerical = zeros(Nx, 1); % Solution vector
    c_numerical(1) = 1; % Initial CO mass fraction fixed at 1.0 at x = 0
    
    % Calculate coefficient for the FTCS scheme
    sigma = D * dt / dx^2;
    
    % Time-stepping loop
    for t = dt:dt:t_final
        % Initialize temporary vector for the new time step
        c_new = c_numerical;
        
        % Loop over space
        for i = 2:Nx - 1
            % Apply FTCS scheme without source term
            c_new(i) = c_numerical(i) + sigma * (c_numerical(i - 1) - 2 * c_numerical(i) + c_numerical(i + 1));
        end
        
        % Apply boundary conditions
        c_new(1) = 1; % CO mass fraction fixed at 1.0 at x = 0
        c_new(end) = 0; % CO mass fraction fixed at 0 at x = L
        
        % Update concentration profile
        c_numerical = c_new;
    end
end