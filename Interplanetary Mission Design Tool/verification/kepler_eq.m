function [output] = kepler_eq(M,E,e)
    %===================================
    % We define the Kepler Equation
    % M is the mean anomaly.
    % e is the eccentricity.
    % E is the aimed solution.
    %===================================
    output = E - e * sin(E) - M;
end