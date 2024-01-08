%=======================================
%
% Interplanetary Mission Analysis
% Design Tool
%
% Master in Space and Aeronautical
% Engineering.
% Astrodynamics Course: Project.
% By: Jorge Simón, Javier Sanchez
% & Iñaki Fernandez.
% Last modification: 20/12/2023
%
%=======================================
clear;
close all;
clc
addpath 'C:\Users\User\Desktop\Iñaki\Fisika\Master\MASE\Courses\Astrodynamics\Scripts\project_code\code_final\assignment1' % Path

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

mu   = 1.32712440018 * 1E20;              % Gravitational parameter [m^3 * s^-2]
Au   = 149597870700;                      % [m]

%=======================================
% Choose departure and arrival planet.
% Mercury = 1; Venus = 2; Earth = 3
% Mars = 4; Jupiter = 5; Saturn = 6
% Uranus = 7; Neptune = 8; Pluto = 9 (extra)
%=======================================

departure_planet = 3;                     % Departure planet index. 
arrival_planet   = 4;                     % Arrival planet index.

%=======================================
% Load orbital elements
%=======================================

planet1_OE = readmatrix(filenames(departure_planet),FileType="text"); % Orbital elements of the departure planet.
planet2_OE = readmatrix(filenames(arrival_planet),  FileType="text"); % Orbital elements of the arrival planet.

%=======================================
% Choose first departure day, departure 
% days quantity and time of flight days
% quantity.
%=======================================

% Departure min day
Y_D                     = 2024;
M_D                     = 1;
D_D                     = 11;
utch_D                  = 0.0;
utcm_D                  = 0.0;
utcs_D                  = 0.0;
[jul_day_D, jul_cent_D] = Julian(Y_D,M_D,D_D,utch_D,utcm_D,utcs_D);

% Departure and TOF day quantity
d_Days                  = 700;
tof_Days                = 700;

%=======================================
% Choose Mode:
%  - Mode = 1 : Ignore inclination (Planar)
%               and eccentricity (Circular)
%  - Mode = 2 : Ignore inclination (Planar)
%  - Mode = 3 : Full mode
%=======================================

% Chose Mode explained above^^
Mode = 3;

%=======================================
% Construc delta V matrix for PCP.
%=======================================

deltaV = zeros(d_Days,tof_Days);

% START LOOPS OVER DEPARTURE DAYS!!!!!!!
for i=1:d_Days

%=======================================
% Define the time arrays for PCP.
%=======================================

[jul_cent_D, jul_cent_A, tof] = Julian_Array(jul_day_D + i-1, tof_Days);

% START LOOPS OVER ARRIVAL DAYS!!!!!!!!!
for j=1:tof_Days

%=======================================
% Compute state vector for the chosen 
% mode.
%=======================================

if Mode == 1
    [r1,v1] = state_vector_JC_PC(planet1_OE,jul_cent_D);
    [r2,v2] = state_vector_JC_PC(planet2_OE,jul_cent_A(j));
elseif Mode == 2
    [r1,v1] = state_vector_JC_P(planet1_OE,jul_cent_D);
    [r2,v2] = state_vector_JC_P(planet2_OE,jul_cent_A(j));
else
    [r1,v1] = state_vector_JC(planet1_OE,jul_cent_D);
    [r2,v2] = state_vector_JC(planet2_OE,jul_cent_A(j));
end

%=======================================
% Lambert Problem Solver.
%
% - tm: Stands for trajectory mode:
%     tm = +1 : (Counter Clock 
%               Wise)
%     tm = -1 : (Clock Wise)
% - n_rev: Stands for revolution quantity
% - zmin: min value of z for z_solver
% - zmax: max value of z for z_solver
%=======================================
tm    = 1;                                 
n_rev = 0;                                 
zmin  = -(pi/2)^2;
zmax  = pi^2;

% Compute the angle
[dtheta] = delta_TrueAnomaly(r1,r2,tm);    % [rads]

% Compute parameters
[A,B,C,P,Q] = params(r1,r2,dtheta);

% Solve transfer function
[z,a,e,p,vt1,vt2] = z_solver_v2(P,Q,tof(j),mu,n_rev, zmin, zmax,dtheta,r1,r2);

%=======================================
% Compute DeltaV.
%=======================================

deltaV1 = vt1 - v1;
deltaV2 = v2 - vt2;
deltaV(i,j) = (abs(norm(deltaV1)) + abs(norm(deltaV2)))*1e-3; % [km/s]

end % END LOOP OF DEPARTURE DAYS!!!!!!!!

loading_bar(i,d_Days,departure_planet,arrival_planet)

end % END LOOP OF ARRIVAL DAYS!!!!!!!!!!

toc

%=======================================
%% POSTPRODUCTION.
%=======================================

PCP_plot(departure_planet, arrival_planet, d_Days, tof_Days, jul_day_D, deltaV, Mode,tof);
