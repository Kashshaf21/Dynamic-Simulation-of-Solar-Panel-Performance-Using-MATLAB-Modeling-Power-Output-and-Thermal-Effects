function eta_t = compute_efficiency(T_vals, eta_ref, beta, T_ref)
% compute_efficiency - Computes the solar panel efficiency over time
%
% Inputs:
%   T_vals   - Vector of cell temperatures over time [°C]
%   eta_ref  - Reference efficiency at T_ref (scalar)
%   beta     - Temperature coefficient (1/°C)
%   T_ref    - Reference temperature [°C]
%
% Output:
%   eta_t    - Efficiency vector over time

    % Efficiency calculation using the given formula
    eta_t = eta_ref * (1 - beta * (T_vals - T_ref));

end
