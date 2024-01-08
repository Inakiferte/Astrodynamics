%=======================================
%
% Cassini's Grand Finale Right Ascen-
% sion and Declination
%
% Master in Space and Aeronautical
% Engineering.
% Astrodynamics Course, Assginment 1.
% By: Jorge Simón Javier Sanchez
% & Iñaki Fernandez.
% Last modification: 29/10/2023
%
%=======================================
clear;
close all;
clc

%=======================================
% Physical Parameters.
%=======================================

G    = 6.673E-11;                         % m^3 / (kg s^2)
mSun = 1.98E30;                           % Sun mass in kg
Au   = 149597870700;                      % m
e    = deg2rad(23.5);                     % Earth inclination

%=======================================
% Casini input data.
% Ref: https://science.nasa.gov/mission/cassini/grand-finale/grand-finale-orbit-guide/
%=======================================

Y    = 2017;
M    = 9;
D    = 15;
utch = 10.0;
utcm = 32.0;
utcs = 0.0;

%=======================================
% Earth orbital elements.
%=======================================

a_E         = 1.00000261;                 % AU
a_dot_E     = 0.00000562;                 % AU / Cy

e_E         = 0.01671123;                 % %
e_dot_E     = -0.00004391;                % 1 / Cy

i_E         = -0.00001531;                % Grad
i_dot_E     = -0.01294668;                % Grad / Cy

Omega_E     = 0.0;                        % Grad
Omega_dot_E = 0.0;                        % Grad / Cy

omega_E     = 102.93768193;               % Grad
omega_dot_E = 0.32327364;                 % Grad / Cy

L_E         = 100.46457166;               % Grad
L_dot_E     = 35999.37244981;             % Grad / Cy

%=======================================
% Saturn orbital elements.
%=======================================

a_S         = 9.53667594;                 % AU
a_dot_S     = -0.00125060;                % AU / Cy

e_S         = 0.05386179;                 % %
e_dot_S     = -0.00050991;                % 1 / Cy

i_S         = 2.48599187;                 % Grad
i_dot_S     = 0.00193609;                 % Grad / Cy

Omega_S     = 113.66242448;               % Grad
Omega_dot_S = -0.28867794;                % Grad / Cy

omega_S     = 92.59887831;                % Grad
omega_dot_S = -0.41897216;                % Grad / Cy

L_S         = 49.95424423;                % Grad
L_dot_S     = 1222.49362201;              % Grad / Cy

%=======================================
% 1/2- Compute Julian Day/Century.
%=======================================

[jul_day, jul_cent] = Julian(Y,M,D,utch,utcm,utcs);

%=======================================
% 3- Propagate the six orbital elements.
%=======================================

a_prop_E     = Propagation(a_E,a_dot_E,jul_cent) * Au;                        % a in meters.
a_prop_S     = Propagation(a_S,a_dot_S,jul_cent) * Au;                        % a in meters.

e_prop_E     = Propagation(e_E,e_dot_E,jul_cent);                             % Earth eccentricity.
e_prop_S     = Propagation(e_S,e_dot_S,jul_cent);                             % Earth eccentricity.

i_prop_E     = mod(deg2rad(Propagation(i_E,i_dot_E,jul_cent)),pi);            % Inclination 0<=i<pi (results in radians).
i_prop_S     = mod(deg2rad(Propagation(i_S,i_dot_S,jul_cent)),pi);            % Inclination 0<=i<pi (results in radians).

Omega_prop_E = mod(deg2rad(Propagation(Omega_E,Omega_dot_E,jul_cent)), 2*pi); % RAAN 0<=Omega<2pi (results in radians).
Omega_prop_S = mod(deg2rad(Propagation(Omega_S,Omega_dot_S,jul_cent)), 2*pi); % RAAN 0<=Omega<2pi (results in radians).

omega_prop_E = mod(deg2rad(Propagation(omega_E,omega_dot_E,jul_cent)), 2*pi); % AP 0<=omega<2pi (results in radians).
omega_prop_S = mod(deg2rad(Propagation(omega_S,omega_dot_S,jul_cent)), 2*pi); % AP 0<=omega<2pi (results in radians).

L_prop_E     = mod(deg2rad(Propagation(L_E,L_dot_E,jul_cent)), 2*pi);         % ML 0<=Omega<2pi (results in radians).
L_prop_S     = mod(deg2rad(Propagation(L_S,L_dot_S,jul_cent)), 2*pi);         % ML 0<=Omega<2pi (results in radians).

%=======================================
% 4- Compute angular momentum.
%=======================================

h_E = Angular_Momentum(G * mSun,a_prop_E,e_prop_E);  % Earth angular momentum.
h_S = Angular_Momentum(G * mSun,a_prop_S,e_prop_S);  % Saturn angular momentum.

%=======================================
% 5- Determine argumnet of perihelion
% and mean anomaly.
%=======================================

omega_E = mod(argp(omega_prop_E,Omega_prop_E), 2*pi);  % Result in Rads
omega_S = mod(argp(omega_prop_S,Omega_prop_S), 2*pi);  % Result in Rads

M_E = mean_anomaly(omega_prop_E, L_prop_E);                  % Result in Rads
M_S = mean_anomaly(omega_prop_S, L_prop_S);                  % Result in Rads


%=======================================
% 6- Solve Kepler Equation
% with Newton-Rhapson method.
%=======================================

N = 100000000;                                   % N-R Max iteration. INTEGER
eps = 10E-12;                                    % N-R accuracy.

E_E = newton_rhapson_method(M_E,e_prop_E,N,eps); % Results in Rads.
E_S = newton_rhapson_method(M_S,e_prop_S,N,eps); % Results in Rads.

%=======================================
% 7- Compute the true anomaly.
%=======================================

theta_E = true_anomaly(e_prop_E,E_E); % Results in Rads.
theta_S = true_anomaly(e_prop_S,E_S); % Results in Rads.

%=======================================
% 8- Orbital elements to state vector.
%=======================================

% Compute geocentric position and velocity values
r_E = zeros(1,3);
v_E = zeros(1,3);

r_S = zeros(1,3);
v_S = zeros(1,3);

cr_E = (h_E^2 / (G * mSun)) * (1.0 / (1.0 + e_prop_E * cos(theta_E)));
cr_S = (h_S^2 / (G * mSun)) * (1.0 / (1.0 + e_prop_S * cos(theta_S)));

r_E(1) = cr_E * cos(theta_E); % + R_S (optional)
r_E(2) = cr_E * sin(theta_E); % + R_S (optional)
r_E(3) = 0.0;                 % + R_S (optional)

r_S(1) = cr_S * cos(theta_S);
r_S(2) = cr_S * sin(theta_S);
r_S(3) = 0.0;

v_E(1) = -(G * mSun / h_E) * sin(theta_E);
v_E(2) = (G * mSun / h_E) * (e_prop_E + cos(theta_E));
v_E(3) = 0.0;

v_S(1) = -(G * mSun / h_S) * sin(theta_S);
v_S(2) = (G * mSun / h_S) * (e_prop_S + cos(theta_S));
v_S(3) = 0.0;

% Compute the rotational matrix
R_S = rotational_matrix(Omega_prop_S, omega_S, i_prop_S);
R_E = rotational_matrix(Omega_prop_E, omega_E, i_prop_E);

% Compute the rotated vectors

r_E_geo = R_E * transpose(r_E);
r_S_geo = R_S * transpose(r_S);
disp(r_S_geo)


% Compute the difference between Saturn and earth
r_geo = r_S_geo - r_E_geo;
r_mod = sqrt(r_geo(1)^2 + r_geo(2)^2 + r_geo(3)^2);

% Compute the celestial latitud and longitud
cel_latitude = asin(r_geo(3) / r_mod);                    % Beta in Rads.
cel_longitud = atan2(r_geo(2), r_geo(1));                 % Lamda in Rads.

% Compute right ascension and declination
declination = asin(sin(cel_latitude) * cos(e) + cos(cel_latitude) * sin(cel_longitud) * sin(e));     % Delta in Rads.
ascension = acos((cos(cel_latitude) * cos(cel_longitud)) / cos(declination));                        % Alpha in Rads.

% Right results in the needed format
[h,m,s] = degrees_to_time(rad2deg(2*pi-ascension));
[degrees_decli, minutes_decli, seconds_decli] = radians_to_dms(declination);

output_file = "results.txt";
file = fopen(output_file, 'w');

fprintf(file, "=======================================================\n");
fprintf(file, "Cassini's Grand Finale Right Ascension and Declination.\n");
fprintf(file, "Master in Space and Aeronautical Engineering.\n");
fprintf(file, "Astrodynamics Course, Assignment 1.\n");
fprintf(file, "By: Jorge Simón, Javier Sanchez & Iñaki Fernandez.\n");
fprintf(file, "=======================================================\n");
fprintf(file, "\n");
fprintf(file, "=======================================================\n");
fprintf(file, "Results for the right ascension and declination:\n");
fprintf(file, " Right ascension alpha: %d h %d min %f s.\n", h, m, s);
fprintf(file, " Right declination delta: %d° %d' %f''.\n", degrees_decli, minutes_decli, seconds_decli);
fprintf(file, "=======================================================\n");

fprintf("=================================================\n");
fprintf("The program has finished successfully!\n");
fprintf("You can find the results in: %s\n", output_file);
fprintf("By: Jorge Simón, Javier Sanchez & Iñaki Fernandez.\n");
fprintf("=================================================\n");

fclose(file);
