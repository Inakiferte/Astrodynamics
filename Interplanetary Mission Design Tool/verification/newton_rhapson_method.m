function [output] = newton_rhapson_method(M,e,N,eps)
    %===================================
    % M is the mean anomaly.
    % e is the esccentricity.
    % N is the iteration steps
    % for N-R method.
    % eps is the accuracy of the method
    %
    % The algorithm:
    %
    % x1 = x0 - f(x0) / f'(x0)
    % Is |x1-x0| < eps?
    % YES: root x1
    % NO : x0 = x1 and repeat.
    %
    % Our functions:
    % f(E) = E - e * sin(E) - M
    % f'(E) = 1 - e * cos(E)
    %===================================

    % Initial guess
    x0 = M;                               % Good initial guess.

    % Iterations for the N-R method
    for i=1: N

        x1 = x0 - (kepler_eq(M,x0,e)) / (kepler_eq_dev(M,x0,e));

        % if for accuracy
        if abs(x1 - x0) < eps
            output = x1;
            break
        else
            x0 = x1;
        end
    end
end