clear;
close all;
clc
tic
%=======================================
% Files containing orbital elements for
% each planet.
%=======================================

% File names
filenames = ["orbital_elements_planets\orbit_elem_Mercury.txt", "orbital_elements_planets\orbit_elem_Venus.txt","orbital_elements_planets\orbit_elem_Earth.txt", ...
    "orbital_elements_planets\orbit_elem_Mars.txt", "orbital_elements_planets\orbit_elem_Jupiter.txt","orbital_elements_planets\orbit_elem_Saturn.txt", ...
    "orbital_elements_planets\orbit_elem_Uranus.txt", "orbital_elements_planets\orbit_elem_Neptune.txt","orbital_elements_planets\orbit_elem_Pluto.txt"];


%=======================================
% Physical Parameters.
%=======================================

G    = 6.673E-11;                         % [m^3 / (kg s^2)]
mSun = 1.98E30;                           % Sun mass in [kg]
mu      = 1.32712440018 * 10^20;          % [m^3 * s^-2]
Au   = 149597870700;                      % [m]

%=======================================
% Choose departure and arrival planet.
% Mercury = 1; Venus = 2; Earth = 3
% Mars = 4; Jupiter = 5; Saturn = 6
% Uranus = 7; Neptune = 8; Pluto = 9 (extra)
%=======================================

departure_planet = 3;                                                 % Departure planet index. It will always be Earth.
arrival_planet   = 4;                                                 % Arrival planet index.

planet1_OE = readmatrix(filenames(departure_planet),FileType="text"); % Orbital elements of the departure planet.
planet2_OE = readmatrix(filenames(arrival_planet),FileType="text");   % Orbital elements of the arrival planet.

%=======================================
% Choose departure and arrival times.
%=======================================

% Departure min day
Y_D                     = 2020;
M_D                     = 1;
D_D                     = 1;
utch_D                  = 0.0;
utcm_D                  = 0.0;
utcs_D                  = 0.0;
[jul_day_D, jul_cent_D] = Julian(Y_D,M_D,D_D,utch_D,utcm_D,utcs_D);


% Departure max day
Y_D_max                 = 2025;
M_D_max                 = 06;
D_D_max                 = 23;
utch_D_max              = 0.0;
utcm_D_max              = 0.0;
utcs_D_max              = 0.0;
[jul_day_D_max, jul_cent_D_max] = Julian(Y_D_max,M_D_max,D_D_max,utch_D_max,utcm_D_max,utcs_D_max);
d_Days                  = jul_day_D_max - jul_day_D;
tof_Days                = 10000;


%=======================================
% Compute the state vector.
% Mode = 1 : Ignore inclination (Planar)
% and eccentricity (Circular)
% Mode = 2 : Ignore inclination (Planar)
% Mode = 3 : Full mode
%=======================================

% Define DeltaV
deltaV = zeros(d_Days,tof_Days);

% Chose Mode explained above^^
Mode = 3;

% START LOOPS OVER DEPARTURE DAYS!!!!!!!
for i=1:d_Days
%=======================================
% Define the time arrays for PCP
%=======================================

[jul_cent_D, jul_cent_A, tof] = Julian_Array(jul_day_D + i - 1, tof_Days);

% START LOOPS OVER ARRIVAL DAYS!!!!!!!!!
for j=1:tof_Days

% Select day and tof for orbits
if i == 200 && j == 269

% % Compute state vector
if Mode == 1
    [r1,v1] = state_vector_JC_PC(planet1_OE,jul_cent_D);
    [r2,v2] = state_vector_JC_PC(planet2_OE,jul_cent_A(j));
elseif Mode == 2
    [r1,v1,R_E,r_E,v_E] = state_vector_JC_P(planet1_OE,jul_cent_D);
    [r2,v2,R2_E,r2_E,v2_E] = state_vector_JC_P(planet2_OE,jul_cent_A(j));
else
    [r1,v1,theta1,a1,e1,Omega1,omega1,i1,r1_E] = state_vector_JC(planet1_OE,jul_cent_D);
    [r2,v2,theta2,a2,e2,Omega2,omega2,i2,r2_E] = state_vector_JC(planet2_OE,jul_cent_A(j));
end

%=======================================
% Obtain Terminal Velocities
% - tm: Stads for trajectory mode:
%     tm = +1 : Short way (Counter Clock 
%               Wise)
%     tm = -1 : Long way (Clock Wise)
% - n_rev: Stands for revolution quantity
% - zmin: min value of z for z_solver
% - zmax: max value of z for z_solver
%=======================================
tm    = 1;                                 % Choose mode
n_rev = 0;                                 % Choose revolution quantity
zmin  = -(pi/2)^2;
zmax  = pi^2;

% Compute the angle
[dtheta] = delta_TrueAnomaly(r1,r2,tm); % [rads]

% Compute Simo solver parameters
[A,B,C,P,Q] = params(r1,r2,dtheta);

% Solve transfer function
 guess = (dtheta / 2)^2;
% guess = 0;
 N = 100000000;                                   % N-R Max iteration. INTEGER
 eps = 10E-12;                                    % N-R accuracy.
%[z,a,vt1,vt2,e,p,ite] = newton_rhapson_method_v2(guess,P,Q,mu,n_rev,tof(j),N,eps,dtheta,r1,r2);
%if ite == 10
 [z,a,e,theta1,Omega,omega,inc,p,vt1,vt2] = z_solver_v2(P,Q,tof(j),mu,n_rev, zmin, zmax,dtheta,r1,r2);
%end

%=======================================
% Compute DeltaV
%=======================================
deltaV1 = vt1 - v1;
deltaV2 = v2 - vt2;
deltaV(i,j) = (abs(norm(deltaV1)) + abs(norm(deltaV2)))*1e-3; %[km/s] 

end % END SELECT DEPARTURE AND TOF

end % END LOOP OF ARRIVAL DAYS!!!!!!!!!!
end % END LOOP OF DEPARTURE DAYS!!!!!!!!

plot_orbit(e1,e2,a1,a2,i1,i2,Omega1,Omega2,omega1,omega2,theta1,theta2,dtheta,a,e,Omega,omega,inc,r1,r2);






