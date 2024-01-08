function [theta] = true_anomaly(e,E)
    % Compute true anomaly.
    theta = 2.0 * atan(tan(E / 2.0) * sqrt((1 + e) / (1 - e)));
end