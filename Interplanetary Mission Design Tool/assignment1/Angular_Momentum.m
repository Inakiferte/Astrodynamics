function [h] = Angular_Momentum(mu,a,e)
    % Compute angular momentum.
    h = sqrt(mu * a * (1 - e^2));
end