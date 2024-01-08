function [R] = rotational_matrix(Omega,omega,inclination)
    % Compute the rotation matrix
    R = zeros(3,3);

    R(1,1) = cos(Omega) * cos(omega) - sin(Omega) * cos(inclination) * sin(omega);
    R(1,2) = -cos(Omega) * sin(omega) - sin(Omega) * cos(inclination) * cos(omega);    
    R(1,3) = sin(Omega) * sin(inclination);
    R(2,1) = sin(Omega) * cos(omega) + cos(Omega) * cos(inclination) * sin(omega);
    R(2,2) = -sin(omega) * sin(Omega) + cos(Omega) * cos(inclination) * cos(omega);    
    R(2,3) = -cos(Omega) * sin(inclination);
    R(3,1) = sin(inclination) * sin(omega);
    R(3,2) = sin(inclination) * cos(omega);
    R(3,3) = cos(inclination);
end