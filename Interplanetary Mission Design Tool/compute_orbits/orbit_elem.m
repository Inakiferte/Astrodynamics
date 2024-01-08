function [a,e1,theta1,Omega,omega,i] = orbit_elem(r1,v1,mu)
%=======================================
% Function that given r1, v1 and mu 
% obtains the orbital elements of the 
% computed trajectory.
%=======================================

% Angular momentum. Don't care if use r2 an v2, is the same
h1 = cross(r1,v1);

% Eccentricity vector and eccentricity
e1_vector = (cross(v1,h1) / mu) - (r1 / norm(r1));
e1 = norm(e1_vector);

% Nodes
N = cross([0 0 1]',h1);

% True anomaly at departure
if dot(r1,v1) >= 0
    theta1 = acos(dot(e1_vector,r1)/(e1 * norm(r1)));
else
    theta1 = 2*pi - acos(dot(e1_vector,r1)/(e1 * norm(r1)));
end

% Inclination
i = acos(h1(3) / norm(h1));

% Longitude of the ascending node
if N(2) >= 0
    Omega = acos(N(1) / norm(N));
else
    Omega = 2*pi - acos(N(1) / norm(N));
end

% Argument of the periapsis
if e1_vector(3) >= 0
    omega = acos(dot(N,e1_vector) / (norm(N) * e1));
else
    omega = 2*pi - acos(dot(N,e1_vector) / (norm(N) * e1));
end

% Semi-major axis
a = 1 / ((2/norm(r1)) - ((norm(v1))^2)/mu);

end