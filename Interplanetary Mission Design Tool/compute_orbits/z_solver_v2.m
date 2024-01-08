function [z,a,e,theta1,Omega,omega,i,p,v_1,v_2] = z_solver_v2(P,Q,tof,mu,n_rev, zmin, zmax,dtetha,r1,r2)
% Function that given r1, r2, tof, dtetha y mu obtains (if it is possible)
% a value for z with a limited minimun error.


% We declare the minimun and maximun 'transfer time'.
dtmin = dt_function(P,Q,zmin,mu,n_rev);
dtmax = dt_function(P,Q,zmax,mu,n_rev);


% It is possible that for the previously defined zmin and zmax values,
% their respective 'transfer time' gives a non-real value and errors appear in the iterative process.
% For this purpose, this value is checked and if it is not real, the value is changed in
% increments of +- 0.1 until a valid result is obtained. For
% zmin we INCREASE its value to move to the right of the graph.
% And for zmax we REDUCE its value to move to the left of the graph.

%a = 0;
while ~isreal(dtmin)
%    %a = a+1;
    zmin = zmin + 0.1;
    dtmin = dt_function(P,Q,zmin,mu,n_rev);
end

%a = 0;

while ~isreal(dtmax)
%    %a = a+1;
    zmax = zmax - 0.1;
    dtmax = dt_function(P,Q,zmax,mu,n_rev); 
end

% Once the zmin and zmax values have been established, the iterative process to find the final z-value is
% using the Bolzano method.
% We calculate zp, which is the mean value between zmin and zmax, to obtain its determined dtp.
% If this dtp is less than tof, zmin becomes zp and if
% dtp is greater than tof, zmax becomes zp. The process is repeated until the minimum accepted difference is obtained.
% The solver is limited to 1000 iterations to avoid entering infinite loops.

comp           = 1;       % Difference between proposed dt and calculated dt.
diferencia_min = 1e-6;    % Minimun accepted difference.
num_vueltas    = 0;       % Iteration counter started at 0.
max_vueltas    = 1000;    % Number of maximun iterations.
t = 1;
while any(comp > diferencia_min) && all(t == 1)

    num_vueltas = num_vueltas + 1;
    zp          = 0.5*(zmin+zmax);
    dtp         = dt_function(P,Q,zp,mu,n_rev);
    comp        = abs(dtp-tof);

    if (dtp < tof)
        zmin = zp;
    else
        zmax = zp;
    end

    if num_vueltas == max_vueltas
        t = 0;
    end
end
z = zp;

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
