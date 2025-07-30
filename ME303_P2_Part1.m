clear; clc; close all;

% Physical properties of the egg and boiling water
k = 0.5; 
rho = 1000; 
c = 3000;
alpha = k/(rho*c); 

% Cooking environment settings
T_s = 100;        
T_target = 80;    
hold_req = 10.0;  

% Numerical grid and time step setup for the heat equation
dr = 5e-4;                    
dt = 0.4*dr^2/alpha;          
max_time = 7200;              

% Egg sizes and initial temperatures

% Egg radii: quail, chicken, ostrich
R_list = [0.012, 0.020, 0.060];    
names   = {'Quail', 'Chicken', 'Ostrich'};
T0 = 5;     

% To store cook times for summary
results = []; 

for e = 1:numel(R_list)
    R = R_list(e);
    r = 0:dr:R;
    N = numel(r);

    % Set initial temperature distribution in the egg
    U = r * T0;
    U(1)   = 0;            % Symmetry at center (Neumann BC)
    U(end) = R * T_s;      % Boundary condition at surface (Dirichlet BC)

    % Prepare variables for time tracking and state logging
    time = 0;
    over = 0;
    cook_time = NaN;
    times = time;
    Tcent = U(2)/dr;       % Estimate center temp

    % For surface plot: store all profiles
    T_all = U ./ r; T_all(1) = U(2)/dr;

  
    while time < max_time
        Unew = U;
        % Heat diffusion update using finite difference
        Unew(2:N-1) = U(2:N-1) + alpha*dt*( U(3:N) - 2*U(2:N-1) + U(1:N-2) )/dr^2;

        % Maintain boundary conditions at each time step
        Unew(1)   = 0;
        Unew(end) = R*T_s;
        U = Unew;
        time = time + dt;

        % Calculate temperature at the center and inside the egg
        Tcenter = U(2)/dr;
        T_interior = U(2:end)./r(2:end);
        Tmin = min([Tcenter, T_interior(1:end-1)]);

        % Record temperature and time data
        times(end+1) = time;
        Tcent(end+1) = Tcenter;

        % Store temperature profiles for visualization
        T_profile = U ./ r; T_profile(1) = Tcenter;
        T_all = [T_all; T_profile];

        % Check if egg temperature is held long enough for doneness
        if Tmin >= T_target
            over = over + dt;
            if over >= hold_req
                cook_time = time;
                break;
            end
        else
            over = 0;
        end
    end

    % Save final cooking time for each egg
    results = [results; e, R, T0, cook_time];

    % Plot center temperature evolution over time
    figure;
    plot(times, Tcent, 'LineWidth', 1); grid on;
    yline(T_target, '--', '80°C target');
    xlabel('Time (s)');
    ylabel('Center temperature (°C)');
    title(sprintf('%s (R=%.3f m), T_0=%g°C', names{e}, R, T0));

    %Visualize temperature across egg radius and time in 3D
    [t_grid, r_grid] = meshgrid(times, r);
    T_all_plot = T_all'; % Transpose for surf so [radius x time]

    figure;
    surf(r_grid, t_grid, T_all_plot, 'EdgeColor', 'none');
    xlabel('Radius (m)');
    ylabel('Time (s)');
    zlabel('Temperature (°C)');
    title(sprintf('3D Surface: %s (R=%.3f m), T_0=%g°C', names{e}, R, T0));
    colorbar; axis tight;
    clim([0, 100]);
    colormap(jet);
    view(45, 30); shading interp;
    drawnow;
end

% Display summary of cooking times for different egg sizes
fprintf('\nCook times (>=%.0f s hold):\n', hold_req);
fprintf('Egg       R (m)   T0 (°C)   Cook time (s)\n');
for k=1:size(results,1)
    eIdx = results(k,1);
    R    = results(k,2);
    T0   = results(k,3);
    tsec = results(k,4);
    fprintf('%-8s  %6.3f    %6.1f     %8.2f\n', names{eIdx}, R, T0, tsec);
end
