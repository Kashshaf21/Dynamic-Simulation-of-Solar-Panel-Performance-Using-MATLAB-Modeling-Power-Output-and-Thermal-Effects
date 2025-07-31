%% PV Power Output Eqyation
function P = pv_power_output(t, eta, G, A)
% Compute PV power output over a time vector using η(t) and G(t) as functions
%
% Inputs:
%   t         : Time vector
%   eta  : Function handle for efficiency η(t)
%   G   : Function handle for irradiance G(t)
%   A         : Panel area [m²]
%
% Output:
%   P         : Power output vector [W]

%eta_vals = eta_func(t);   % Evaluate efficiency at each time point
%G_vals = G_func(t);       % Evaluate irradiance at each time point
P = eta .* A .* G;
end