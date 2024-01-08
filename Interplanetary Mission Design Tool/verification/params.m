function [A,B,C,P,Q] = params(r1,r2,delta_theta)
r1_n = norm(r1);
r2_n = norm(r2);

A = sqrt(r1_n);
B = sqrt(r2_n) * cos(delta_theta/2);
C = sqrt(r2_n) * sin(delta_theta/2);
P = A^2 + B^2+ C^2;
Q = 2 * A * B;

end