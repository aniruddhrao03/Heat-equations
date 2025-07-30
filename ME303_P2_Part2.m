clear; clc; close all;

% Problem parameters
alpha   = 2;                   
Nx      = 50;      

%Grid Settings
dx      = 1/Nx;                % Spatial grid size
dt      = 0.5*dx^2/alpha;      % Explicit time step (stability)
x       = linspace(0,1,Nx+1);  % Spatial grid

% Fourier series terms for analytical solution
Nterms  = 100;                  
times   = [0.001, 0.01, 0.1, 10]; % Times for comparison

% Grid independence study parameters
Nxs       = [25, 50, 100, 200, 400]; % Test grid sizes
t_GI      = 0.1;                     % Time for error comparison
CFLfac    = 0.25;                    % CFL scaling factor

% Initialize Vectors for Errors
dx_vec    = zeros(size(Nxs));
err_inf   = zeros(size(Nxs));
err_L2    = zeros(size(Nxs));

% Numerical time-integration setup 
max_t      = max(times);             % Maximum time of interest
Nsteps_max = round(max_t/dt);        % Number of steps to reach max_t
U_all      = zeros(Nx+1, Nsteps_max+1); % Store all time steps for surface plots
U_all(:,1) = cos(pi*x);              % Initial condition
timevec    = zeros(1, Nsteps_max+1); 

% Explicit time-stepping for main solution
for step = 1:Nsteps_max
    u_old = U_all(:,step);
    u_new = u_old;
    for i = 2:Nx
        u_new(i) = u_old(i) + alpha*dt/dx^2 * (u_old(i+1) - 2*u_old(i) + u_old(i-1));
    end

    % Dirichlet boundary conditions
    u_new(1) = 0; 
    u_new(end) = 2;
    U_all(:,step+1) = u_new;
    timevec(step+1) = step*dt;
end

% Analytical vs numerical comparisons at selected times 
for k = 1:length(times)
    t  = times(k);
    xa = linspace(0,1,100);     % Fine grid for plotting
    ua = 2*xa;                  % Analytical steady-state

    for n = 1:Nterms            % Analytical transient terms
        if mod(n,2)==1
            bn = -4/(n*pi);
        else
            bn = (2*((1/(n-1)) + (1/(n+1))) + 4/n)/pi;
        end
        ua = ua + bn*exp(-alpha*n^2*pi^2*t) .* sin(n*pi*xa);
    end
    idx = round(t/dt);          % Match time step index

    figure;
    plot(xa, ua, '-b', 'LineWidth', 1.4); 
    hold on;

    plot(x, U_all(:, idx+1), 'r', 'LineWidth', 1);

    legend('Analytical', 'Numerical');
    title(['Comparison at t = ', num2str(t)]);
    xlabel('x'); ylabel('u(x,t)'); grid on;   
end

%Grid independence study: solve for various dx and compute errors 
for k = 1:numel(Nxs)
    Nxk = Nxs(k); dxk = 1/Nxk;
    dtk = CFLfac * dxk^2 / alpha;
    xk  = linspace(0,1,Nxk+1).';
    nsteps = ceil(t_GI/dtk);
    u = cos(pi*xk);     % Initial condition

    for step = 1:nsteps
        u_old = u; u_new = u_old;
        for i = 2:Nxk
            u_new(i) = u_old(i) + alpha*dtk/dxk^2 * (u_old(i+1) - 2*u_old(i) + u_old(i-1));
        end
        u_new(1) = 0; u_new(end) = 2;
        u = u_new;
    end


    % Analytical reference on same grid
    u_exact = 2*xk;
    for n = 1:Nterms
        if mod(n,2)==1
            bn = -4/(n*pi);
        else
            bn = (2*((1/(n-1)) + (1/(n+1))) + 4/n)/pi;
        end
        u_exact = u_exact + bn * exp(-alpha*n^2*pi^2*t_GI) .* sin(n*pi*xk);
    end
    dx_vec(k)  = dxk;
    err_inf(k) = norm(u - u_exact, inf);           % Max norm error
    err_L2(k)  = sqrt(trapz(xk, (u - u_exact).^2)); % L2 error
end

%Error convergence plots (log–log scale) 

figure;
loglog(dx_vec, err_inf, 'b', 'LineWidth', 1.3); grid on;
xlabel('\Delta x'); ylabel('||error||_{\infty}');
title('Grid Independence: L_\infty error vs \Delta x');

figure;
loglog(dx_vec, err_L2, 'r', 'LineWidth', 1.3); grid on;
xlabel('\Delta x'); ylabel('L_2 error');
title('Grid Independence: L_2 error vs \Delta x');


