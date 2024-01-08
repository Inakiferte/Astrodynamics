function [jul_cent_D, jul_cent_A, tof] = Julian_Array(jul_day_D, tof_Days)
%=======================================
% This function computes the required
% time arrays for the Lambert Problem.
%
% INPUTS:
%       - jul_day_D : Departure julian
%                     day.
%       - tof_Days  : Time of flight
%                     days.
% OUTPUTS: 
%       - jul_cent_D : Departure day
%                      julian century.
%       - jul_cent_A : Array with arrival
%                      days for a given
%                      departure day.
%       - tof        : Array of time of 
%                      flights for a 
%                      given departure.
%=======================================

% Define arrays
tof = zeros(1,tof_Days);
jul_cent_A = zeros(1,tof_Days);

% Compute tof days
for i=1:tof_Days
    tof(i) = i * 86400;
end

% Julian centuries
jul_cent_D = (jul_day_D - 2451545) / 36525;

for i=1:tof_Days
    jul_cent_A(i) = ((jul_day_D + i) - 2451545) / 36525;
end

end