function [a,e,p] = orbit_elem(z,r1,r2,delta_theta)
[A,B,C,P,Q] = params(r1, r2, delta_theta);
a = (P - Q * c_function(0, z))/(2 * z * c_function(1,z)^2);
e = sqrt(1 - (A^2 * C^2 / (a^2 * z * c_function(1,z)^2))) ;
p = (2 * A^2 * C^2) / (P - (Q * c_function(0,z)));
end