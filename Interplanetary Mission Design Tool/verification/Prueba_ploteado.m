clc
clear all

%% Ploteado delta t vs z

%% Datos y ploteado grafica pg 48

r1 = [0 1 0]; % Random value of R1
r2 = [1 0 0]; % Random value of R2
mu = 1; % Gravitational parameter SUN
tof = 0:0.2:20;
dtetha1 = pi/2;
dtetha2 = 2*pi/3;
dtetha3 = 3*pi/2;
dtetha4 = 5*pi/3;
z1 = zeros(1, length(tof));
z2 = zeros(1, length(tof));
z3 = zeros(1, length(tof));
z4 = zeros(1, length(tof));
zmin = -(pi/2)^2;
zmax = pi^2;
r1n = norm(r1);
r2n = norm(r2);
P = r1n + r2n;
Q1 = 2 * sqrt(r1n * r2n) * cos(dtetha1 / 2);
Q2 = 2 * sqrt(r1n * r2n) * cos(dtetha2 / 2);
Q3 = 2 * sqrt(r1n * r2n) * cos(dtetha3 / 2);
Q4 = 2 * sqrt(r1n * r2n) * cos(dtetha4 / 2);


for i=1:length(tof)

    [z1(i)] = z_solver(P,Q1,tof(i),mu,0, zmin, zmax,dtetha1,r1,r2);
    [z2(i)] = z_solver(P,Q2,tof(i),mu,0, zmin, zmax,dtetha2,r1,r2);
    [z3(i)] = z_solver(P,Q3,tof(i),mu,0, zmin, zmax,dtetha3,r1,r2);
    [z4(i)] = z_solver(P,Q4,tof(i),mu,0, zmin, zmax,dtetha4,r1,r2);

end

% Suponiendo que z1, z2, z3 y z4 son conjuntos de datos que se están graficando
% Encuentra los valores máximos de los conjuntos z1, z2, z3 y z4
max_z_values = [max(z1), max(z2), max(z3), max(z4)];

% Encuentra el máximo de todos los valores máximos anteriores
max_z = max(max_z_values);

% Graficar con líneas más gruesas y límite en el eje x desde -2 hasta el máximo de los datos
plot(z1, tof, 'b', 'LineWidth', 2); % Graficar y1 en azul con línea gruesa
hold on; % Mantener la gráfica actual
plot(z2, tof, 'r', 'LineWidth', 2); % Graficar y2 en rojo con línea gruesa
plot(z3, tof, 'y', 'LineWidth', 2); % Graficar y3 en amarillo con línea gruesa
plot(z4, tof, 'm', 'LineWidth', 2); % Graficar y4 en magenta con línea gruesa

xlim([-2, max_z]); % Establecer límite en el eje x desde -2 hasta el máximo de los datos

grid on;

xlim([-(pi/2)^2, (pi)^2]);

xticks([-(pi/2)^2, 0, (pi/2)^2, (pi)^2]);
xticklabels({'-(\pi/2)^2', '0', '(\pi/2)^2', '\pi^2'});

% Añadir líneas verticales en los puntos específicos (-(pi/2)^2, 0, (pi/2)^2 y pi^2)
xline(-(pi/2)^2, '--k', 'LineWidth', 1); % Línea vertical en x = -(pi/2)^2
xline(0, '--k', 'LineWidth', 1); % Línea vertical en x = 0
xline((pi/2)^2, '--k', 'LineWidth', 1); % Línea vertical en x = (pi/2)^2
xline((pi)^2, '--k', 'LineWidth', 1); % Línea vertical en x = pi^2

xlabel('Z (NU)');
ylabel('\Delta t (NU)');
title('\Delta t vs Z plot');
legend('\Delta\theta = pi/2', '\Delta\theta = 2*pi/3', '\Delta\theta = 3*pi/2', '\Delta\theta = 5*pi/3'); % Leyenda de las funciones

hold off; 

