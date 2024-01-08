clear;
close all;
clc

% Define same input as D. De la Torre et al.
% All in NDU.
r1      = [1,0,0];
r2      = [0,-2,0];
delta_t = linspace (0, 50, 200);
z       = zeros(1, length(delta_t));
mu      = 1;

% Define revolution value
nrev = 1;

% Define delta_theta values
delta_theta = [nrev * 2*pi + 3*pi/2];

if nrev>0
    zmin = 0.00001;
    zmax = pi^2;
    for j=1:length(delta_theta)
        [A,B,C,P,Q] = params(r1,r2,delta_theta(j));
    for i=1:length(delta_t)
        [z(i),a,e,p,v_1,v_2] = z_solver(P,Q,delta_t(i),mu,nrev,zmin,zmax,delta_theta(j),r1,r2);
    end


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
    plot(z_space,dt_direct, 'b', 'LineWidth', 4);
    xlabel('Z (NU)');
ylabel('\Delta t (NU)');
title('\Delta t vs Z plot');
    hold on
    axis([-(pi/2)^2 pi^2  0 50]);
    x_lower = -(pi/2)^2;
x_upper = (pi)^2;

x_lower_limit = -(pi/2)^2;
x_upper_limit = pi^2;

x_mid_point = (pi^2 + x_lower_limit) / 2 + 0.08 * (pi^2 - x_lower_limit); 


x_ticks = [x_lower_limit, 0, x_mid_point, x_upper_limit];


axis([x_lower_limit, x_upper_limit, 0, 50]);


x_tick_labels = {'-π^2', '0', '(π/2)^2', 'π^2'}; 
xticks(x_ticks);
xticklabels(x_tick_labels);


grid_values = [x_lower_limit, 0, x_mid_point, x_upper_limit];
grid on;
set(gca, 'xtick', grid_values);

y_val = 8 * pi;
yline(y_val, '--k', 'LineWidth', 1); 


text(x_lower + 0.01* (x_upper - x_lower), y_val + 2, '\Delta t = 8\pi', 'FontSize', 15); 
text(x_lower + 0.1 * (x_upper - x_lower), y_val-20 , '\Delta \theta = 2\pi + \pi/2', 'FontSize', 15);

    end
else
    zmin = -(pi/2)^2;
    zmax = pi^2;
    for j=1:length(delta_theta)
        [A,B,C,P,Q] = params(r1,r2,delta_theta(j));
    for i=1:length(delta_t)
        [z(i),a,e,p,v_1,v_2] = z_solver(P,Q,delta_t(i),mu,nrev,zmin,zmax,delta_theta(j),r1,r2);
    end


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
    plot(z_space,dt_direct, 'b', 'LineWidth', 4);
    xlabel('Z (NU)');
ylabel('\Delta t (NU)');
title('\Delta t vs Z plot');
    hold on
    axis([-(pi/2)^2 + 0.01 pi^2  0 50]);

    end
end



