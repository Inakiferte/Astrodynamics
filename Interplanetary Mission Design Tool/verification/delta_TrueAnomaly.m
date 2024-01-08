function [dtheta] = delta_TrueAnomaly(r1,r2,tm)
%=======================================
% It computes the angle between r1 and r2
% tm: Stads for trajectory mode
%   tm = +1 : Short way (Counter Clock Wise)
%   tm = -1 : Long way (Clock Wise)
%=======================================

% Angle with the cosine
dtheta_p   = acos((dot(r1,r2)) / (norm(r1) * norm(r2)));

% Ensure right cuadrant for the chosen mode
dtheta     = asin(tm * sqrt(1 - cos(dtheta_p)^2));       % Output in [rads]

end