function plot_orbit(e1,e2,a1,a2,i1,i2,Omega1,Omega2,omega1,omega2,theta1,theta2,dtheta,a,e,Omega,omega,inc,r1_dep,r2_arr)


AU = 1.496*1e11;   %[m]
mu = 1.32712440018 * 10^20;  % mu - gravitational parameter of central body

%=======================================
% Orbital elements.
%=======================================

h1_E = Angular_Momentum(mu,a1,e1);  % Planet angular momentum.
h2_E = Angular_Momentum(mu,a2,e2);  % Planet angular momentum.

% Compute position and velocity values
r1_E = zeros(1,3);
v1_E = zeros(1,3);
r2_E = zeros(1,3);
v2_E = zeros(1,3);

% Theta Planets
theta_planet_1 = linspace (theta1, theta1 + 2*pi, 5000);
theta_planet_2 = linspace (theta2, theta2 + 2*pi, 5000);

% Compute the rotational matrix
R1_E = rotational_matrix(Omega1, omega1, i1);
R2_E = rotational_matrix(Omega2, omega2, i2);

r1 = zeros(length(theta_planet_1),3);
r2 = zeros(length(theta_planet_2),3);
for i=1:length(theta_planet_1)
cr1_E(i) = (h1_E^2 / (mu)) * (1.0 / (1.0 + e1 * cos(theta_planet_1(i))));
cr2_E(i) = (h2_E^2 / (mu)) * (1.0 / (1.0 + e2 * cos(theta_planet_2(i))));

r1_E(i,1) = cr1_E(i) * cos(theta_planet_1(i)); 
r1_E(i,2) = cr1_E(i) * sin(theta_planet_1(i)); 
r1_E(i,3) = 0.0; 

% Compute the rotated vectors
r1(i,:) = R1_E * (r1_E(i,:))';


r2_E(i,1) = cr2_E(i) * cos(theta_planet_2(i)); 
r2_E(i,2) = cr2_E(i) * sin(theta_planet_2(i)); 
r2_E(i,3) = 0.0; 

% Compute the rotated vectors
r2(i,:) = R2_E * (r2_E(i,:))';

end


% Trajectory
h_E = Angular_Momentum(mu,a,e);  % Planet angular momentum.
theta_trajectory = linspace (theta1, theta1 + dtheta, 5000);

R_E = rotational_matrix(Omega, omega, inc);

r = zeros(length(theta_trajectory),3);
r_E = zeros(length(theta_trajectory),3);

for i=1:(length(theta_trajectory))
cr_E(i) = (h_E^2 / (mu)) * (1.0 / (1.0 + e * cos(theta_trajectory(i))));

r_E(i,1) = cr_E(i) * cos(theta_trajectory(i)); 
r_E(i,2) = cr_E(i) * sin(theta_trajectory(i)); 
r_E(i,3) = 0.0; 

% Compute the rotated vectors
r(i,:) = R_E * (r_E(i,:))';

end






% Plot
plot3(r1(:,1)/AU,r1(:,2)/AU,r1(:,3)/AU,'LineWidth', 2, 'Color','b');
hold on
plot3(r2(:,1)/AU,r2(:,2)/AU,r2(:,3)/AU,'LineWidth', 2, 'Color','r');
plot3(r(:,1)/AU,r(:,2)/AU,r(:,3)/AU,'LineWidth', 2, 'Color','g');
plot3(r1_dep(1)/AU,r1_dep(2)/AU,r1_dep(3)/AU,'o', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
plot3(r2_arr(1)/AU,r2_arr(2)/AU,r2_arr(3)/AU,'o', 'MarkerFaceColor', 'r', 'MarkerSize', 15);

zlim([-1.5, 1.5]);
xlabel('X [AU]','Interpreter','latex',FontSize=25)
ylabel('Y [AU]','Interpreter','latex',FontSize=25)
zlabel('Z [AU]','Interpreter','latex',FontSize=25)
leyenda = legend('Earth', 'Mars', 'Spacecraft');
set(gcf, 'color', 'k'); % 'k' for black background, 'none' for no color
set(gca, 'color', 'k');
set(gca, 'xcolor', 'white');
set(gca, 'ycolor', 'white');
set(gca, 'zcolor', 'white');
set(leyenda, 'Color', 'none','Interpreter','latex', 'TextColor', 'white', 'Location', 'North West'); 
set(gca, 'FontSize', 25); 


end