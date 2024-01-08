function [dtheta] = delta_TrueAnomaly(r1,r2,tm)
%=======================================
% It computes the difference in true 
% anomaly between r1 and r2.
%
% INPUTS:
%       - tm: Stads for trajectory mode
%           tm = +1 : (Counter Clock Wise)
%           tm = -1 :  (Clock Wise)
%       - r1: Position array of planet 1
%       - r2: Position array of planet 2
% OUTPUT:
%       - dtheta: Delta true anomaly in
%                 in rads!!!!!!!!!!!!!!.
%=======================================

% Find cosine
ca = (dot(r1,r2)) / (norm(r1) * norm(r2));

% Find sine
sa = sqrt(1 - (ca)^2);

% Compute dtheta with the tangent
% This ensures the right cuadrant!
dtheta = atan2(sa, ca);

% Compute the z direction
vect = cross(r1,r2);

% If z is negative, impose Counter
% Clock Wise
if vect(3) < 0
    dtheta = 2*pi - dtheta;
end

% Change direction if we choose Clock Wise
if tm == -1
    dtheta = 2*pi-dtheta;
end

end