%=======================================
%
% Z solver verification.
%
% Master in Space and Aeronautical
% Engineering.
% Astrodynamics Course: Project.
% By: Jorge Simón, Javier Sanchez
% & Iñaki Fernandez.
% Last modification: 03/12/2023
%
%=======================================
clear;
close all;
clc

% Define same input as D. De la Torre et al.
% All in NDU.
r1      = [1,0,0];
r2      = [0,-1,0];
delta_t = linspace (0, 20, 200);
z       = zeros(1, length(delta_t));
mu      = 1;

% Define revolution value
nrev = 1;

% Define delta_theta values
delta_theta = [nrev * 2*pi + pi/2,nrev * 2*pi + 2*pi/3,nrev * 2*pi + 3*pi/2,nrev * 2*pi + 5*pi/3];

if nrev>0
    zmin = 0.00001;
    zmax = pi^2;
    for j=1:length(delta_theta)
        [A,B,C,P,Q] = params(r1,r2,delta_theta(j));
    for i=1:length(delta_t)
        [z(i),a,e,p,v_1,v_2] = z_solver(P,Q,delta_t(i),mu,nrev,zmin,zmax,delta_theta(j),r1,r2);
    end

    figure(1);
    plot(z,delta_t, '-', 'LineWidth', 1);
    axis([-(pi/2)^2 pi^2  0 20]);
    legend('$2\pi + \pi /2$', '$2\pi + 2\pi / 3$', '$2\pi + 3\pi / 2$', '$2\pi + 5\pi / 3$', 'Interpreter', 'latex','FontSize', 15,'Location', 'SouthEast');
    xlabel('$z$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('$\Delta t$', 'Interpreter', 'latex', 'FontSize', 15);
    hold on
    end
    % Direct deltat
    z_space  = linspace(0.00001,pi^2,2000);
    dt_direct = zeros(1, length(z_space));
    for j=1:length(delta_theta)
    [A,B,C,P,Q] = params(r1,r2,delta_theta(j));
    for i=1: length(z_space)
        dt = dt_function(P,Q,z_space(i),mu,nrev);
        dt_direct(i) = dt;
    end
    
    figure(2);
    plot(z_space,dt_direct, '-', 'LineWidth', 1);
    xlabel('$z$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('$\Delta t$ Direct', 'Interpreter', 'latex', 'FontSize', 15);
    legend('$2\pi + \pi /2$', '$2\pi + 2\pi / 3$', '$2\pi + 3\pi / 2$', '$2\pi + 5\pi / 3$', 'Interpreter', 'latex','FontSize', 15,'Location', 'SouthEast');
    hold on
    axis([-(pi/2)^2 pi^2  0 20]);
    end
else
    zmin = -(pi/2)^2;
    zmax = pi^2;
    for j=1:length(delta_theta)
        [A,B,C,P,Q] = params(r1,r2,delta_theta(j));
    for i=1:length(delta_t)
        [z(i),a,e,p,v_1,v_2] = z_solver(P,Q,delta_t(i),mu,nrev,zmin,zmax,delta_theta(j),r1,r2);
    end

    figure(1);
    plot(z,delta_t, '-', 'LineWidth', 1);
    xlim([-(pi/2)^2 + 0.01, pi^2 ]);
    legend('$\pi /2$', '$2\pi / 3$', '$3\pi / 2$', '$5\pi / 3$', 'Interpreter', 'latex','FontSize', 15);
    xlabel('$\textbf{z}$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('\textbf{$\Delta t$}', 'Interpreter', 'latex', 'FontSize', 15);
    hold on
    end
    % Direct deltat
    z_space  = linspace(-2,8.5,2000);
    dt_direct = zeros(1, length(z_space));
    for j=1:length(delta_theta)
    [A,B,C,P,Q] = params(r1,r2,delta_theta(j));
    for i=1: length(z_space)
        dt = dt_function(P,Q,z_space(i),mu,nrev);
        dt_direct(i) = dt;
    end
    
    figure(2);
    plot(z_space,dt_direct, '-', 'LineWidth', 1);
    xlabel('$z$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('$\Delta t$ Direct', 'Interpreter', 'latex', 'FontSize', 15);
    legend('$\pi /2$', '$2\pi / 3$', '$3\pi / 2$', '$5\pi / 3$', 'Interpreter', 'latex','FontSize', 15);
    hold on
    axis([-(pi/2)^2 + 0.01 pi^2  0 20]);
    end
end

%=======================================
%% Derivative Stuff
%=======================================
clc
% Plot of the derivative
if nrev>0
z_space  = linspace(0.00001,pi^2,2000);
dt_deriv = zeros(1, length(z_space));
x = [0,4];
y = [0,0];
for j=1:length(delta_theta)
for i=1: length(z_space)
    [A,B,C,P,Q] = params(r1,r2,delta_theta(j));
    [dt_deriv(i)] = dt_derivative_v2(P,Q,z_space(i),mu,nrev);
end
figure(3);
plot(z_space,dt_deriv, '-', 'LineWidth', 1);
hold on;
xlabel('$z$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\Delta t/dz$', 'Interpreter', 'latex', 'FontSize', 15);
legend('$2\pi + \pi /2$', '$2\pi + 2\pi / 3$', '$2\pi + 3\pi / 2$', '$2\pi + 5\pi / 3$', 'Interpreter', 'latex','FontSize', 15,'Location', 'SouthEast');
axis([0 4  -25 10]);
end
plot(x,y,'--',LineWidth=1,Color='green');

% Find roots
[z_root] = zero_finder(delta_theta,nrev,mu,r1,r2);

% Display the result
for i=1: length(delta_theta)
    fprintf('For a variation of true anomaly of %f\n', delta_theta(i))
    fprintf('The zero of the derivative is at z = %f\n', z_root(i));
end
end

%=======================================
%% NR Solver tests for n_rev > 0
%=======================================
clc
for j = 1:length(delta_theta)
[A,B,C,P,Q] = params(r1,r2,delta_theta(j));
N = 100000000;                                   % N-R Max iteration. INTEGER
eps = 10E-15;                                    % N-R accuracy.
min = z_root(j);
tof = linspace(dt_function(P,Q,min,mu,nrev), 20, 200); 
z_first = zeros(1,length(tof));
z_second = zeros(1,length(tof));
% Time of Flight, delta t

% Solve
for i = 1:length(tof)
[output] = newton_rhapson_method_nrev(min,P,Q,mu,nrev,tof(i),N,eps);
z_first(i) = output(1);
z_second(i) = output(2);
disp(output)
end
figure(4);
plot(z_first,tof, '-', 'LineWidth', 1);
hold on;
xlabel('$z$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\Delta t/dz$', 'Interpreter', 'latex', 'FontSize', 15);
legend('$2\pi + \pi /2$', '$2\pi + 2\pi / 3$', '$2\pi + 3\pi / 2$', '$2\pi + 5\pi / 3$', 'Interpreter', 'latex','FontSize', 15,'Location', 'SouthEast');
axis([-(pi/2)^2 pi^2  0 tof(end)]);
plot(z_second,tof, '-', 'LineWidth', 1);
end

[a_test,e_test,p_test] = orbit_elem(0.45,r1,r2,delta_theta(4));
d_theta_test = linspace (0, 2*pi, 500);
for i=1:length(d_theta_test)
r(i) = (a_test * (1-e_test^2))/(1 + e_test * cos(d_theta_test(i)));
x(i) = r(i) * cos(d_theta_test(i));
y(i) = r(i) * sin(d_theta_test(i));
end
[a_test_2,e_test_2,p_test_2] = orbit_elem(0.98,r1,r2,delta_theta(4));
%d_theta_test = linspace (0, 2*pi, 500);
for i=1:length(d_theta_test)
r_2(i) = (a_test_2 * (1-e_test_2^2))/(1 + e_test_2 * cos(d_theta_test(i)));
x_2(i) = r_2(i) * cos(d_theta_test(i));
y_2(i) = r_2(i) * sin(d_theta_test(i));
end
figure(5)
plot(x,y,'o-', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold on
plot(x_2,y_2,'o-','MarkerFaceColor', 'r', 'MarkerSize', 5);
