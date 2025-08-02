clear; clc; close all;

%Properties
k = 0.64;  
rho = 1000;  
c = 4180;
alpha = k/(rho*c);

%Cooking settings
T_surface = 100;    
T_target  = 80;      
hold_req  = 10;   

% Grid / time
dr = 5e-4;  
 
dt = 0.4*dr^2/alpha;

%For While Loop logic
max_time = 7200;               

% Eggs and Initial Conditions
R_list = [0.012, 0.020, 0.060];                
names  = {'Quail','Chicken','Ostrich'};
T_init = 5;                                    
results = zeros(numel(R_list),4);              % [idx, R, T0, cook_time]

for e = 1:numel(R_list)
    R = R_list(e);
    r = (0:dr:R).';  N = numel(r);

    % IC
    T = T_init*ones(N,1);
    T(end) = T_surface;

    % Logs
    time = 0.0; over_duration = 0.0; cook_time = NaN;
    times = time; Tcent = T(1);
    T_all = T.';  % rows = time snapshots

    %Time loop
    while time < max_time
        T = forward_step(T, alpha, r, dr, dt, T_surface);
        time = time + dt;

        % Logs
        times(end+1,1) = time;        
        Tcent(end+1,1) = T(1);        
        T_all(end+1,:) = T.';         

        
        % Check if the minimum temperature is above the target temperature
        if min(T(1:end-1)) >= T_target
            % Accumulate the duration above the target temperature
            over_duration = over_duration + dt;
            % Check if the required hold time has been reached
            if over_duration >= hold_req
                cook_time = time; break
            end
        else
            % Reset the over duration if the target temperature is not met
            over_duration = 0.0;
        end
    end

    results(e,:) = [e, R, T_init, cook_time];

    % Plot center temperature
    figure('Color','w');
    plot(times, Tcent, 'LineWidth', 1.4); grid on;
    yline(T_target,'--','80°C target');
    xlabel('Time (s)'); ylabel('Center temperature (°C)');
    title(sprintf('%s (R=%.3f m), T_0=%g°C', names{e}, R, T_init));

    % 3D profile
    [tt, rr] = meshgrid(times, r);
    Z = T_all.';  % radius × time
    figure('Color','w');
    surf(rr, tt, Z, 'EdgeColor','none');
    xlabel('Radius (m)'); ylabel('Time (s)'); zlabel('Temperature (°C)');
    title(sprintf('T(r,t) — %s (R=%.3f m)', names{e}, R));
    colorbar; clim([0 100]); colormap(jet);
    view(45,30); shading interp; drawnow;
end

%Summary
fprintf('\nCook times (≥ %.0f°C for ≥ %.0f s):\n', T_target, hold_req);
fprintf('Egg       R (m)   T0 (°C)   Cook time (s)   Cook time (min)\n');
for k = 1:size(results,1)
    idx = results(k,1); Rk = results(k,2); T0k = results(k,3); tk = results(k,4);
    if isnan(tk)
        fprintf('%-8s  %6.3f    %6.1f     %12s   %14s\n', names{idx}, Rk, T0k, '—', '—');
    else
        fprintf('%-8s  %6.3f    %6.1f     %12.2f   %14.2f\n', names{idx}, Rk, T0k, tk, tk/60);
    end
end


% FTCS step on T (spherical) 
function T_next = forward_step(T, alpha, r, dr, dt, T_surface)
% One explicit step: dT/dt = α( T_rr + (2/r)T_r ), central differences
    N = numel(r);
    T_next = T;
    for i = 2:(N-1)
        d2 = (T(i+1)-2*T(i)+T(i-1))/dr^2;
        d1 = (T(i+1)-T(i-1))/(2*dr);
        T_next(i) = T(i) + dt*alpha*( d2 + (2/r(i))*d1 );
    end
    T_next(1)   = T_next(2);     
    T_next(end) = T_surface;     
end
