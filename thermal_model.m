function [T,time]=thermal_model(C_th, R_th, alpha, T0, time_array, G, Tamb)
    % thermal_model Simulates thermal behavior of a PV panel
    % Parameters:
    % C_th  - Thermal capacitance [J/°C]
    % R_th  - Thermal resistance [°C/W]
    % alpha - Fraction of irradiance converted to heat
    % T0    - Initial panel temperature [°C]
    % time_array - Time array [s]
    % G     - Array of irradiance values G(t) [W/m²]
    % Tamb  - Array of ambient temperature values Tamb(t) [°C]

    % Use the provided time array
    time = time_array;
    N = length(time);
    T = zeros(1, N);
    T(1) = T0;
    
    % Ensure G and Tamb are the same length as time
    if length(G) ~= N
        error('G array must have the same length as time array. G length: %d, time length: %d', length(G), N);
    end
    if length(Tamb) ~= N
        error('Tamb array must have the same length as time array. Tamb length: %d, time length: %d', length(Tamb), N);
    end

    % Calculate time step
    dt = time(2) - time(1);

    % Simulation loop using Euler's method
    for i = 1:N-1
        G_t = G(i);
        Tamb_t = Tamb(i);
        dTdt = (1/C_th) * (alpha * G_t - (T(i) - Tamb_t) / R_th);
        T(i+1) = T(i) + dt * dTdt;
    end

    % Plot PV panel temperature
    %figure;
    %plot(time/3600, T, 'b', 'LineWidth', 2);
    %xlabel('Time [hours]');
    %ylabel('Panel Temperature [°C]');
    %title('PV Panel Temperature Over Time');
    %grid on;
end