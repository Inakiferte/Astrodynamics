function [E_prop] = Propagation(E,E_dot,jul_cent)
    % Propagation of orbital elements.
    E_prop = E + E_dot * jul_cent;
end