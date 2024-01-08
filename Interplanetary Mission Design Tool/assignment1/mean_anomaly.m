function [M] = mean_anomaly(omega_prop,L_prop)
    % Compute mean anomaly
    M = L_prop - omega_prop;
end