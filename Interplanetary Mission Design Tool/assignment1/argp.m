function [omega] = argp(omega_prop,Omega_prop)
    % Compute argument perihelion and mean anomaly.
    omega = omega_prop - Omega_prop;
end