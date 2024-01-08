function [output] = dt_function_NR(P,Q,z,mu,n_rev,dt)
% Function that given the values of P,Q,z y mu computes the 'trasnfer time' (dt)
% It is defined in small pieces named as 'aux' to avoid typos
R = 2 * P * c_function(3,4 * z) + Q * (c_function(1,z) * c_function(2,4 * z) - 2 * c_function(0,z) * c_function(3,4 * z));
S = (c_function(1,z))^3;
V = P - Q * c_function(0,z);
T = sqrt(2 * V / mu);
U = (pi * n_rev/(c_function(1,z))^3) * sqrt((1/(2 * mu)) * ((P - Q * c_function(0,z)) / z)^3);



output = (R/S) * T + U - dt;
end