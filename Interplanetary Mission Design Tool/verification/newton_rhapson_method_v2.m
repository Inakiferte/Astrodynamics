function [output,a,v_1,v_2,e,p] = newton_rhapson_method_v2(guess,P,Q,mu,n_rev,dt,N,eps,dtetha,r1,r2)
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
    % Our functions: (For example Kepler's equation
    % f(E) = E - e * sin(E) - M
    % f'(E) = 1 - e * cos(E)
    %===================================

    % Initial guess
    x0 = guess;                               % Good initial guess.

    % Iterations for the N-R method
    for i=1: N

        x1 = x0 - (dt_function_NR(P,Q,x0,mu,n_rev,dt)) / (dt_derivative_NR(P,Q,x0,mu,n_rev));

        % if for accuracy
        if abs(x1 - x0) < eps
            output = x1;
            break
        else
            x0 = x1;
        end
    end
z = x1;

[a,e,p] = orbit_elem(z,r1,r2,dtetha);

r1_n = norm(r1);
r2_n = norm(r2);

f = 1 - (r1_n/p)*(1-cos(dtetha));
g_d = f;
g = ((r1_n*r2_n)/sqrt(mu*p))*sin(dtetha);
    
v_1 = (r2-f*r1)/g;
v_2 = (g_d*r2 - r1)/g;
end