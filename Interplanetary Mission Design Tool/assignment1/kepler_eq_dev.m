function [output] = kepler_eq_dev(M,E,e)
    %===================================
    % We define the Kepler Equation
    % derivative.
    % E is the aimed solution.
    %===================================
    output = 1.0 - e * cos(E);
end