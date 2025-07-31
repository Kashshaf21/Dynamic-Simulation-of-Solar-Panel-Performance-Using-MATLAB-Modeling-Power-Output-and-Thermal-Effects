%% Energy Output Integration

function E = energy_output(P, t, tol)
% Computes energy output using Adaptive Simpson's Rule
%
% Inputs:
%   P      : Array of power values P(t_i)
%   t      : Array of time points corresponding to P values
%   tol    : Desired accuracy (e.g. 1e-4)
%
% Output:
%   E      : Total energy output (e.g. in Joules if time in seconds)

    if nargin < 3
        tol = 1e-8;
    end
    
    % Input validation
    if length(P) ~= length(t)
        error('P and t arrays must have the same length');
    end
    
    if length(P) < 2
        error('At least 2 data points are required');
    end
    
    % Sort by time if not already sorted
    [t, idx] = sort(t);
    P = P(idx);
    
    % Create interpolation function
    P_func = @(time) interp1(t, P, time, 'linear', 'extrap');
    
    % Get time interval bounds
    a = t(1);
    b = t(end);
    
    maxDepth = 20; % prevent infinite recursion
    E = recurse(P_func, a, b, tol, simpson(P_func, a, b), maxDepth);

end

% Simpson's rule on a single interval
function S = simpson(f, a, b)
    c = (a + b) / 2;
    S = (b - a) / 6 * (f(a) + 4 * f(c) + f(b));
end

% Recursive adaptive Simpson's method
function S = recurse(f, a, b, tol, whole, depth)
    c = (a + b) / 2;
    left = simpson(f, a, c);
    right = simpson(f, c, b);
    if depth <= 0 || abs(left + right - whole) < 15 * tol
        S = left + right + (left + right - whole) / 15;  % correction term
    else
        S = recurse(f, a, c, tol/2, left, depth - 1) + ...
            recurse(f, c, b, tol/2, right, depth - 1);
    end
end
