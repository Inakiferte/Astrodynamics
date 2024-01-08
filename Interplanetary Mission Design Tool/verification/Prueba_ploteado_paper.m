clc
clear all

%% Ploteado delta t vs z

%% Datos y ploteado grafica pg 48

r1 = [1 0 0]; % Random value of R1
r2 = [-2 0 0]; % Random value of R2
mu = 1; % Gravitational parameter SUN
tof = 0:0.2:30;
dtetha1 = pi/2;
dtetha2 = 3*pi/2;
% dtetha3 = 3*pi/2;
% dtetha4 = 5*pi/3;
z1 = zeros(1, length(tof));
z2 = zeros(1, length(tof));
% z3 = zeros(1, length(tof));
% z4 = zeros(1, length(tof));
zmin = -(pi/2)^2;
zmax = pi^2;
r1n = norm(r1);
r2n = norm(r2);
P = r1n + r2n;
Q1 = 2 * sqrt(r1n * r2n) * cos(dtetha1 / 2);
Q2 = 2 * sqrt(r1n * r2n) * cos(dtetha2 / 2);


for i=1:length(tof)

    [z1(i)] = z_solver(P,Q1,tof(i),mu,0, zmin, zmax,dtetha1,r1,r2);
    [z2(i)] = z_solver(P,Q2,tof(i),mu,0, zmin, zmax,dtetha2,r1,r2);

end


mid_value = (pi)^2 / 4; % 
shift_amount = 0.0000001 * ((pi)^2 - 0); 
mid_x_position = ((pi)^2 + 0) / 2 + shift_amount; 

x_lower = -(pi/2)^2;
x_upper = (pi)^2;

new_xticks = [x_lower, 0, mid_x_position, x_upper];

xticklabels_new = {'-(\pi/2)^2', '0', '(\pi/2)^2', '\pi^2'};

% Plot con los cambios
figure (1);
plot(z1, tof, 'b', 'LineWidth', 4); 
hold on;
xlim([x_lower, x_upper]); 
grid on;

xticks(new_xticks); 
xticklabels(xticklabels_new); 


y_val = 2 * pi;
yline(y_val, '--k', 'LineWidth', 1); 


text(x_lower + 0.1 * (x_upper - x_lower), y_val + 2, '\Delta t = 2\pi', 'FontSize', 18); % Ajuste en y
text(x_lower + 0.1 * (x_upper - x_lower), y_val + 20, '\Delta \theta = \pi/2', 'FontSize', 18);


xlabel('Z (NU)');
ylabel('\Delta t (NU)');
title('\Delta t vs Z plot');
hold off;

%% Otro plot

figure(2);
plot(z2, tof, 'b', 'LineWidth', 4);
xlim([x_lower, x_upper]); 
grid on;

xticks(new_xticks); 
xticklabels(xticklabels_new); 


y_val = 2 * pi;
yline(y_val, '--k', 'LineWidth', 1); 


text(x_lower + 0.1 * (x_upper - x_lower), y_val + 2, '\Delta t = 2\pi', 'FontSize', 18);
text(x_lower + 0.1 * (x_upper - x_lower), y_val + 20, '\Delta \theta = 3\pi/2', 'FontSize', 18);


xlabel('Z (NU)');
ylabel('\Delta t (NU)');
title('\Delta t vs Z plot');
hold off;






