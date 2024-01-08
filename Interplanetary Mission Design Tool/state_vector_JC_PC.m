function [r,v] = state_vector_JC_PC(planet_OE,jul_cent)
%=======================================
% From the Orbital Elements of a planet
% and the Julian Century it computes
% the state vector.
% P : Stands for planar orbits, i = 0
% C : Stands for circular orbits, e = 0
%
% INPUTS:
%       - planet_OE: Matrix with the 
%                    orbital elements.
%       - jul_cent : Julian century to
%                    which we want to
%                    propagate the 
%                    orbital elements.
% OUTPUTS:
%       - r        : position vector.
%       - v        : velocity vector.
%=======================================

%=======================================
% Physical Parameters.
%=======================================

mu      = 1.32712440018 * 1E20;          % [m^3 * s^-2]  
Au      = 149597870700;                  % m

%=======================================
% Orbital elements.
%=======================================

a_E         = planet_OE(1,1);            % AU
a_dot_E     = planet_OE(2,1);            % AU / Cy

e_E         = 0.0;                       % #
e_dot_E     = 0.0;                       % 1 / Cy

i_E         = 0.0;                       % Grad
i_dot_E     = 0.0;                       % Grad / Cy

Omega_E     = planet_OE(1,4);            % Grad
Omega_dot_E = planet_OE(2,4);            % Grad / Cy

omega_E     = planet_OE(1,5);            % Grad
omega_dot_E = planet_OE(2,5);            % Grad / Cy

L_E         = planet_OE(1,6);            % Grad
L_dot_E     = planet_OE(2,6);            % Grad / Cy


%=======================================
% 3- Propagate the six orbital elements.
%=======================================

a_prop_E     = Propagation(a_E,a_dot_E,jul_cent) * Au;                        % a in meters.

e_prop_E     = Propagation(e_E,e_dot_E,jul_cent);                             % Earth eccentricity.

i_prop_E     = mod(deg2rad(Propagation(i_E,i_dot_E,jul_cent)),pi);            % Inclination 0<=i<pi (results in radians).

Omega_prop_E = mod(deg2rad(Propagation(Omega_E,Omega_dot_E,jul_cent)), 2*pi); % RAAN 0<=Omega<2pi (results in radians).

omega_prop_E = mod(deg2rad(Propagation(omega_E,omega_dot_E,jul_cent)), 2*pi); % AP 0<=omega<2pi (results in radians).

L_prop_E     = mod(deg2rad(Propagation(L_E,L_dot_E,jul_cent)), 2*pi);         % ML 0<=Omega<2pi (results in radians).

%=======================================
% 4- Compute angular momentum.
%=======================================

h_E = Angular_Momentum(mu,a_prop_E,e_prop_E);  % Planet angular momentum.

%=======================================
% 5- Determine argumnet of perihelion
% and mean anomaly.
%=======================================

omega_E = mod(argp(omega_prop_E,Omega_prop_E), 2*pi);  % Result in Rads

M_E = mean_anomaly(omega_prop_E, L_prop_E);            % Result in Rads

%=======================================
% 6- Solve Kepler Equation
% with Newton-Rhapson method.
%=======================================

N = 100000000;                                   % N-R Max iteration. INTEGER
eps = 10E-12;                                    % N-R accuracy.

E_E = newton_rhapson_method(M_E,e_prop_E,N,eps); % Results in Rads.

%=======================================
% 7- Compute the true anomaly.
%=======================================

theta_E = true_anomaly(e_prop_E,E_E);            % Results in Rads.

%=======================================
% 8- Orbital elements to state vector.
%=======================================

% Compute position and velocity values
r_E = zeros(1,3);
v_E = zeros(1,3);

cr_E = (h_E^2 / (mu)) * (1.0 / (1.0 + e_prop_E * cos(theta_E)));

r_E(1) = cr_E * cos(theta_E); 
r_E(2) = cr_E * sin(theta_E); 
r_E(3) = 0.0;                 

v_E(1) = -(mu / h_E) * sin(theta_E);
v_E(2) = (mu / h_E) * (e_prop_E + cos(theta_E));
v_E(3) = 0.0;

% Compute the rotational matrix
R_E = rotational_matrix(Omega_prop_E, omega_E, i_prop_E);

% Compute the rotated vectors
r = R_E * transpose(r_E);
v = R_E * transpose(v_E);
end