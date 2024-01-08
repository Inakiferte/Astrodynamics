function [z,a,e,theta1,Omega,omega,i,p,v_1,v_2,ite] = newton_rhapson_method_v2(guess,P,Q,mu,n_rev,dt,N,eps,dtetha,r1,r2)
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
    for ite=1: N
            if ite == 10
                z = 0;
                a = 0;
                v_1 = 0;
                v_2 = 0;
                e = 0;
                p = 0;
                theta1 = 0;
                Omega = 0;
                omega = 0;
                i = 0;
                break
            end
            x1 = x0 - (dt_function_NR(P,Q,x0,mu,n_rev,dt)) / (dt_derivative_NR(P,Q,x0,mu,n_rev));

            % if for accuracy
            if abs(x1 - x0) < eps
                z = x1;
                break
            else
                x0 = x1;
            end
           
        %end
    end
if ite < 10
z = x1;

[A,B,C,P,Q] = params(r1, r2, dtetha);

r1_n = norm(r1);
r2_n = norm(r2);

p = (2 * A^2 * C^2) / (P - (Q * c_function(0,z)));

% Lagrange Parameters
f = 1 - (r2_n/p)*(1-cos(dtetha));
g_d = 1 - (r1_n / p) * (1 - cos(dtetha));
g = ((r1_n*r2_n)/sqrt(mu*p))*sin(dtetha);

% Terminal Velocities
v_1 = (r2-f*r1)/g;
v_2 = (g_d*r2 - r1)/g;

[a,e,theta1,Omega,omega,i] = orbit_elem(r1,v_1,mu);

end
end