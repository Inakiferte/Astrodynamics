function [output] = dt_derivative_NR(P,Q,z,mu,n_rev)
% From D. De La Torre et al.: On the solution of Lambert's problem by regularization
R = 2 * P * c_function(3,4 * z) + Q * (c_function(1,z) * c_function(2,4 * z) - 2 * c_function(0,z) * c_function(3,4 * z));
S = (c_function(1,z))^3;
V = P - Q * c_function(0,z);
T = sqrt(2 * V / mu);
U = (pi * n_rev/(c_function(1,z))^3) * sqrt((1/(2 * mu)) * ((P - Q * c_function(0,z)) / z)^3);

% De La Torre (89)-(97)
V = P - (Q * c_function(0,z));
dc0_dz = (-c_function(1,z)) / 2;
dc1_dz = (c_function(3,z) - c_function(2,z)) / 2;
dc24z_dz = 2 * (2 * c_function(4,4*z) - c_function(3,4*z));
dc34z_dz = 2 * (3 * c_function(5,4*z) - c_function(4,4*z));
dRdz = 2 * P * dc34z_dz + Q * (c_function(2,4*z) * dc1_dz + c_function(1,z) * dc24z_dz - 2 * c_function(3,4*z) * dc0_dz - 2 * c_function(0,z) * dc34z_dz);
dSdz = 3 * c_function(1,z)^2 * dc1_dz;
dTdz = (-Q / sqrt(2 * mu * V)) * dc0_dz;
dUdz = (-3 * n_rev * pi) / (2 * z^2 * c_function(1,z)^4) * sqrt(V / (2 * z * mu)) * (z * c_function(1,z) * Q * dc0_dz + V * c_function(1,z) + 2 * z * V * dc1_dz);
output = (T / S) * dRdz - (R * T / S^2) * dSdz + (R / S) * dTdz + dUdz;
end