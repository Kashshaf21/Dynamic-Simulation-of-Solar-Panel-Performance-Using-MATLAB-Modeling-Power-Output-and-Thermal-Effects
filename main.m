
%% Solar Panel Performance Modeling - Simplified Version (No Cloud/Wind Effects)
% thermal behavior and electrical performance of a PV panel
% using first-order ODE for temperature dynamics
clc;
clear all;
close all;

%% Panel Specifications and Constants

T_ref = 25;             % Reference temperature (°C)
%P_rated           % Rated power at STC (W)
%G_ref             % Reference irradiance (W/m²)
%absorptance        % Absorptance coefficient (dimensionless)
%R_th             % Thermal resistance (K·m²/W)
%C_th            % Thermal capacitance (J/K·m²) 

% Additional panel parameters for I-V characteristics
%V_oc_ref         % Open circuit voltage at STC (V)
%I_sc_ref         % Short circuit current at STC (A)
V_mp_ref = 32.0;        % Maximum power point voltage at STC (V)
I_mp_ref = 9.38;        % Maximum power point current at STC (A)
% Temperature coefficients for I-V characteristics
alpha_I = 0.0005;       % Temperature coefficient of current (A/°C)
beta_V = -0.0032;       % Temperature coefficient of voltage (V/°C)





% Launch the app
app = tutorialApp;

% Pause execution until user clicks Analyze
uiwait(app.UIFigure);



% Show it
%fprintf('Panel Area A = %.2f m²\n', A);



    
%% Environmental Conditions (Time Series)
% Time vector (24 hours with fine resolution for ODE)
t = 0:0.01:24;          % Time in hours (fine resolution for ODE)
dt = t(2) - t(1);       % Time step in hours
dt_sec = dt * 3600;     % Time step in seconds
n_points = length(t);

% Solar irradiance profile (W/m²) 
G0 = 1000 * max(0, sin(pi * (t - 6) / 12)) .* (t >= 6 & t <= 18);


% Season adjustments
        switch season
            case 'Summer'
                season_factor = 1.0;
                temp_base = 25;
            case 'Spring'
                season_factor = 0.8;
                temp_base = 20;
            case 'Winter'
                season_factor = 0.6;
                temp_base = 10;
        end
        
        % Location adjustments
        switch location
            case 'Tropical'
                location_temp_offset = 5;
            case 'Temperate'
                location_temp_offset = 0;
            case 'Cold'
                location_temp_offset = -10;
        end
        
        % Weather adjustments
        switch weather
            case 'Clear Sky'
                weather_factor = 1.0;
            case 'Partly Cloudy'
                weather_factor = 0.7;
            case 'Overcast'
                weather_factor = 0.3;
        end
        
        % Final irradiance
        G = G0 * season_factor * weather_factor;
        
        % Ambient temperature profile
        T_amb = temp_base + location_temp_offset + ...
                8 * sin(pi * (t - 6) / 12) + ...
                3 * sin(2 * pi * (t - 6) / 12);




%% 1. Thermal Modeling Using First-Order ODE
fprintf('=== THERMAL MODELING WITH ODE ===\n');
T0 = T_amb(1);

% Convert time to seconds for the thermal model
t_sec = t * 3600;  % Convert hours to seconds

% Call thermal model with proper time array
[T_cell, time_sec] = thermal_model(C_th, R_th, absorptance, T0, t_sec, G, T_amb);

% Convert time back to hours for consistency with the rest of the code
time = time_sec / 3600;


%% 2. Efficiency Calculation
fprintf('=== EFFICIENCY CALCULATION ===\n');

eta_ref = P_rated / (G_ref * A);
% Calculate efficiency using: eta(t) = eta_ref [1 - beta* (T(t) - T_ref)]
eta =  compute_efficiency(T_cell, eta_ref, beta, T_ref);

% Ensure efficiency doesn't go negative
eta = max(eta, 0);

% Calculate efficiency losses due to temperature
eta_losses = eta_ref - eta;


%% 3. Power Output Calculation
fprintf('=== POWER OUTPUT CALCULATION ===\n');

% Calculate power output using: P(t) = eta(t) × A × G(t)
P_output = pv_power_output(t, eta, G, A);

% Calculate I-V characteristics for comparison
V_oc = zeros(size(t));
I_sc = zeros(size(t));
V_mp = zeros(size(t));
I_mp = zeros(size(t));
P_iv = zeros(size(t));
FF = zeros(size(t));  % Fill factor

for i = 1:n_points
    if G(i) > 0
        % Temperature difference
        dT = T_cell(i) - T_ref;
        
        % Temperature-corrected parameters
        V_oc(i) = V_oc_ref + beta_V * dT;
        I_sc(i) = I_sc_ref * (G(i) / G_ref) * (1 + alpha_I * dT);
        V_mp(i) = V_mp_ref + beta_V * dT;
        I_mp(i) = I_mp_ref * (G(i) / G_ref) * (1 + alpha_I * dT);
        
        % Power at maximum power point
        P_iv(i) = V_mp(i) * I_mp(i);
        
        % Fill factor calculation
        FF(i) = P_iv(i) / (V_oc(i) * I_sc(i));
    end
end

%% 4. Energy Output Integration
fprintf('=== ENERGY OUTPUT INTEGRATION ===\n');

% Daily energy calculation using trapezoidal integration
E_daily = energy_output(P_output, t,1e-8);  % Daily energy in kWh

% Cumulative energy throughout the day
E_cumulative = cumtrapz(t, P_output) / 1000;  % Cumulative energy in kWh


%% 5. Performance Metrics and Analysis
fprintf('=== PERFORMANCE METRICS ===\n');

% Ideal energy without temperature effects
P_ideal = P_rated * G / G_ref;
E_ideal = trapz(t, P_ideal) / 1000;

% Performance ratio
PR = E_daily / E_ideal;

% Temperature losses
temp_losses_energy = trapz(t, eta_losses .* A .* G) / 1000;

% Thermal time constant
tau_thermal = C_th * R_th;

% Maximum temperature gradient
dT_dt_array = [0, diff(T_cell) / dt_sec];
dT_dt_max = max(abs(dT_dt_array));


%% 6. Main Visualization
fprintf('=== GENERATING VISUALIZATIONS ===\n');

% Create comprehensive figure
fig1 = figure('Name', 'Solar Panel Performance Analysis', 'Position', [50, 50, 1200, 800]);

% Subplot 1: Environmental conditions
subplot(2, 4, 1);
yyaxis left
plot(t, G, 'r-', 'LineWidth', 2);
ylabel('Irradiance (W/m²)', 'Color', 'r');
ylim([0, 1200]);
yyaxis right
plot(t, T_amb, 'b-', 'LineWidth', 2);
ylabel('Ambient Temp (°C)', 'Color', 'b');
xlabel('Time (hours)');
title('Environmental Conditions');
grid on;

% Subplot 2: Cell temperature vs ambient
subplot(2, 4, 2);
plot(t, T_cell, 'k-', 'LineWidth', 2);
hold on;
plot(t, T_amb, 'b--', 'LineWidth', 1.5);
ylabel('Temperature (°C)');
xlabel('Time (hours)');
title('Cell vs Ambient Temperature');
legend('Cell Temp', 'Ambient Temp', 'Location', 'best');
grid on;

% Subplot 3: Efficiency variation
subplot(2, 4, 3);
plot(t, eta * 100, 'g-', 'LineWidth', 2);
ylabel('Efficiency (%) / Loss (%)');
xlabel('Time (hours)');
title('Efficiency and Temperature Losses');
legend('Efficiency',  'Location', 'best');
grid on;

% Subplot 4: Power output comparison
subplot(2, 4, 4);
plot(t, P_output, 'r-', 'LineWidth', 2);
hold on;
plot(t, P_iv, 'b--', 'LineWidth', 1.5);
plot(t, P_ideal, 'g:', 'LineWidth', 1.5);
ylabel('Power (W)');
xlabel('Time (hours)');
title('Power Output Comparison');
legend('eta×A×G Method', 'I-V Method', 'Ideal Power', 'Location', 'best');
grid on;

% Subplot 5: Temperature gradient
subplot(2, 4, 5);
plot(t, dT_dt_array, 'k-', 'LineWidth', 1.5);
ylabel('dT/dt (°C/s)');
xlabel('Time (hours)');
title('Temperature Gradient');
grid on;

% Subplot 6: I-V characteristics
subplot(2, 4, 6);
yyaxis left
plot(t, V_oc, 'b-', 'LineWidth', 1.5);
ylabel('Voltage (V)', 'Color', 'b');
yyaxis right
plot(t, I_sc, 'r-', 'LineWidth', 1.5);
ylabel('Current (A)', 'Color', 'r');
xlabel('Time (hours)');
title('I-V Characteristics');
grid on;

% Subplot 7: Cumulative energy
subplot(2, 4, 7);
plot(t, E_cumulative, 'b-', 'LineWidth', 2);
ylabel('Cumulative Energy (kWh)');
xlabel('Time (hours)');
title('Energy Accumulation');
grid on;

% Subplot 8: Fill factor
subplot(2, 4, 8);
plot(t, FF, 'k-', 'LineWidth', 2);
ylabel('Fill Factor');
xlabel('Time (hours)');
title('Fill Factor Variation');
grid on;

axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.5, 0.98, 'Solar Panel Performance Analysis (Simplified Model)', ...
     'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

% Make sure figure 1 is visible
figure(fig1);


%% 8. Parametric Analysis
fprintf('=== GENERATING PARAMETRIC ANALYSIS ===\n');

fig3 = figure('Name', 'Parametric Analysis', 'Position', [150, 150, 1000, 600]);

% Temperature coefficient sensitivity
subplot(2, 3, 1);
beta_range = 0.002:0.0005:0.006;
E_daily_beta = zeros(size(beta_range));
for i = 1:length(beta_range)
    eta_temp = eta_ref * (1 - beta_range(i) * (T_cell - T_ref));
    eta_temp = max(eta_temp, 0);
    P_temp = eta_temp .* A .* G;
    E_daily_beta(i) = trapz(t, P_temp) / 1000;
end
plot(beta_range, E_daily_beta, 'b-', 'LineWidth', 2);
xlabel('Temperature Coefficient ? (1/°C)');
ylabel('Daily Energy (kWh)');
title('Effect of Temperature Coefficient');
grid on;

% Thermal capacitance effect
subplot(2, 3, 2);
C_th_range = 20000:5000:80000;
max_temp_C = zeros(size(C_th_range));
for i = 1:length(C_th_range)
    T_temp = zeros(size(t));
    T_temp(1) = T_amb(1);
    for j = 1:n_points-1
        dT_dt = (1/C_th_range(i)) * (absorptance * G(j) - (T_temp(j) - T_amb(j))/R_th);
        T_temp(j+1) = T_temp(j) + dT_dt * dt_sec;
    end
    max_temp_C(i) = max(T_temp);
end
plot(C_th_range/1000, max_temp_C, 'r-', 'LineWidth', 2);
xlabel('Thermal Capacitance (kJ/K·m²)');
ylabel('Maximum Temperature (°C)');
title('Effect of Thermal Capacitance');
grid on;

% Thermal resistance effect
subplot(2, 3, 3);
R_th_range = 0.02:0.005:0.08;
max_temp_R = zeros(size(R_th_range));
for i = 1:length(R_th_range)
    T_temp = zeros(size(t));
    T_temp(1) = T_amb(1);
    for j = 1:n_points-1
        dT_dt = (1/C_th) * (absorptance * G(j) - (T_temp(j) - T_amb(j))/R_th_range(i));
        T_temp(j+1) = T_temp(j) + dT_dt * dt_sec;
    end
    max_temp_R(i) = max(T_temp);
end
plot(R_th_range, max_temp_R, 'g-', 'LineWidth', 2);
xlabel('Thermal Resistance (K·m²/W)');
ylabel('Maximum Temperature (°C)');
title('Effect of Thermal Resistance');
grid on;

% Absorptance effect
subplot(2, 3, 4);
alpha_range = 0.7:0.025:0.95;
max_temp_alpha = zeros(size(alpha_range));
for i = 1:length(alpha_range)
    T_temp = zeros(size(t));
    T_temp(1) = T_amb(1);
    for j = 1:n_points-1
        dT_dt = (1/C_th) * (alpha_range(i) * G(j) - (T_temp(j) - T_amb(j))/R_th);
        T_temp(j+1) = T_temp(j) + dT_dt * dt_sec;
    end
    max_temp_alpha(i) = max(T_temp);
end
plot(alpha_range, max_temp_alpha, 'k-', 'LineWidth', 2);
xlabel('Absorptance ?');
ylabel('Maximum Temperature (°C)');
title('Effect of Absorptance');
grid on;

% Panel area effect
subplot(2, 3, 5);
A_range = 1:0.2:4;
E_daily_area = zeros(size(A_range));
for i = 1:length(A_range)
    P_temp = eta .* A_range(i) .* G;
    E_daily_area(i) = trapz(t, P_temp) / 1000;
end
plot(A_range, E_daily_area, 'm-', 'LineWidth', 2);
xlabel('Panel Area (m²)');
ylabel('Daily Energy (kWh)');
title('Effect of Panel Area');
grid on;

% Reference efficiency effect
subplot(2, 3, 6);
eta_ref_range = 0.10:0.01:0.25;
E_daily_eta_ref = zeros(size(eta_ref_range));
for i = 1:length(eta_ref_range)
    eta_temp = eta_ref_range(i) * (1 - beta * (T_cell - T_ref));
    eta_temp = max(eta_temp, 0);
    P_temp = eta_temp .* A .* G;
    E_daily_eta_ref(i) = trapz(t, P_temp) / 1000;
end
plot(eta_ref_range * 100, E_daily_eta_ref, 'c-', 'LineWidth', 2);
xlabel('Reference Efficiency (%)');
ylabel('Daily Energy (kWh)');
title('Effect of Reference Efficiency');
grid on;


axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.5, 0.98, 'Parametric Analysis of Key Parameters', ...
     'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

% Make sure figure 3 is visible
figure(fig3);

%% 9. Monthly Analysis
fprintf('=== GENERATING MONTHLY ANALYSIS ===\n');

fig4 = figure('Name', 'Monthly Analysis', 'Position', [200, 200, 1000, 600]);

% Monthly solar radiation data (kWh/m²/day)
monthly_radiation = [2.5, 3.5, 4.8, 6.2, 7.1, 7.5, 7.8, 7.2, 5.8, 4.2, 2.8, 2.2];
monthly_avg_temp = [5, 8, 12, 16, 20, 24, 26, 25, 21, 16, 10, 6];  % °C
months = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
 
% Calculate monthly energy production
monthly_energy = zeros(size(monthly_radiation));
for i = 1:12
    % Average efficiency for the month
    eta_monthly = eta_ref * (1 - beta * (monthly_avg_temp(i) - T_ref));
    eta_monthly = max(eta_monthly, 0);
    
    % Monthly energy (assuming 30 days average)
    monthly_energy(i) = eta_monthly * A * monthly_radiation(i) * 30;
end

annual_energy = sum(monthly_energy);

% Plot monthly analysis
subplot(2, 2, 1);
bar(1:12, monthly_energy);
set(gca, 'XTick', 1:12, 'XTickLabel', months);
xlabel('Month');
ylabel('Energy (kWh)');
title('Monthly Energy Production');
grid on;

subplot(2, 2, 2);
bar(1:12, monthly_radiation);
set(gca, 'XTick', 1:12, 'XTickLabel', months);
xlabel('Month');
ylabel('Solar Radiation (kWh/m²/day)');
title('Monthly Solar Radiation');
grid on;

subplot(2, 2, 3);
plot(1:12, monthly_avg_temp, 'ro-', 'LineWidth', 2, 'MarkerSize', 6);
set(gca, 'XTick', 1:12, 'XTickLabel', months);
xlabel('Month');
ylabel('Temperature (°C)');
title('Monthly Average Temperature');
grid on;

subplot(2, 2, 4);
monthly_efficiency = zeros(size(monthly_avg_temp));
for i = 1:12
    monthly_efficiency(i) = eta_ref * (1 - beta * (monthly_avg_temp(i) - T_ref));
end
plot(1:12, monthly_efficiency * 100, 'go-', 'LineWidth', 2, 'MarkerSize', 6);
set(gca, 'XTick', 1:12, 'XTickLabel', months);
xlabel('Month');
ylabel('Efficiency (%)');
title('Monthly Average Efficiency');
grid on;

axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.5, 0.98, 'Annual Performance Projection', ...
     'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

% Make sure figure 4 is visible
figure(fig4);

% Display all figures
fprintf('Generated 4 figures:\n');
fprintf('Figure 1: Solar Panel Performance Analysis\n');
fprintf('Figure 2: I-V Characteristics Analysis\n');
fprintf('Figure 3: Parametric Analysis\n');
fprintf('Figure 4: Monthly Analysis\n');

% Optional: pause to view each figure
fprintf('\nPress Enter to continue or Ctrl+C to stop...\n');
pause;
% Performance summary
fprintf('Daily Energy Output: %.2f kWh\n', E_daily);
fprintf('Ideal Energy Output: %.2f kWh\n', E_ideal);
fprintf('Peak Power Output: %.2f W\n', max(P_output));
fprintf('Average Cell Temperature: %.1f °C\n', mean(T_cell(G > 0)));
fprintf('Maximum Cell Temperature: %.1f °C\n', max(T_cell));
fprintf('Temperature Rise at Peak: %.1f °C\n', max(T_cell) - T_ref);
fprintf('Performance Ratio: %.3f\n', PR);
fprintf('Temperature Losses: %.2f kWh (%.1f%%)\n', temp_losses_energy, temp_losses_energy/E_ideal*100);
fprintf('Thermal Time Constant: %.1f seconds\n', tau_thermal);
fprintf('Maximum Temperature Gradient: %.4f °C/s\n', dT_dt_max);
fprintf('Average Efficiency: %.1f%%\n', mean(eta(G > 0)) * 100);
fprintf('Peak Efficiency: %.1f%%\n', max(eta) * 100);



%% 10. Summary Report
fprintf('\n=== PERFORMANCE SUMMARY REPORT ===\n');
fprintf('=====================================\n');
fprintf('THERMAL ANALYSIS:\n');
fprintf('- Maximum Cell Temperature: %.1f °C\n', max(T_cell));
fprintf('- Temperature Rise: %.1f °C\n', max(T_cell) - T_ref);
fprintf('- Thermal Time Constant: %.1f seconds\n', tau_thermal);
fprintf('- Max Temperature Gradient: %.4f °C/s\n', dT_dt_max);

fprintf('\nEFFICIENCY ANALYSIS:\n');
fprintf('- Reference Efficiency: %.1f%%\n', eta_ref * 100);
fprintf('- Average Efficiency: %.1f%%\n', mean(eta(G > 0)) * 100);
fprintf('- Peak Efficiency: %.1f%%\n', max(eta) * 100);
fprintf('- Temperature Coefficient: %.1f%%/°C\n', beta * 100);

fprintf('\nENERGY PERFORMANCE:\n');
fprintf('- Daily Energy Output: %.2f kWh\n', E_daily);
fprintf('- Ideal Energy (no temp effect): %.2f kWh\n', E_ideal);
fprintf('- Performance Ratio: %.3f\n', PR);
fprintf('\nPOWER CHARACTERISTICS:\n');
fprintf('- Peak Power Output: %.2f W\n', max(P_output));
fprintf('- Rated Power: %.2f W\n', P_rated);
fprintf('- Peak Power Factor: %.3f\n', max(P_output)/P_rated);

fprintf('\nANNUAL PROJECTION:\n');
fprintf('- Estimated Annual Energy: %.1f kWh\n', annual_energy);
fprintf('- Capacity Factor: %.1f%%\n', annual_energy/(P_rated*8760/1000)*100);

fprintf('\nSIMULATION DETAILS:\n');
fprintf('- Time Step: %.2f minutes\n', dt * 60);
fprintf('- Total Simulation Points: %d\n', n_points);
fprintf('- Panel Area: %.1f m²\n', A);
fprintf('=====================================\n');

fprintf('\n=== SIMULATION COMPLETE ===\n');
fprintf('Generated 4 comprehensive analysis figures\n');