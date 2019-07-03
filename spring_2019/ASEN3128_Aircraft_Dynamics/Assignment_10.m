% ASEN 3218
% Assignment 9
% Matthew Ryan
% April 8, 2019

% Housekeeping
close all; clear; clc;

% Flight parameters (Boeing 747 Case II)
g = 9.81; % gravitaional acceleration | m/s^2
rho = .652694; % atmospheric density | kg/m^3
W = 6.366*10^5 * 4.44822162; % 747 weight | N
u0 = 518 * .3048; % nominal airspeed | m/s
xi = -6.8*pi/180; % | rad
theta0 = 0; % climb angle | rad
b = 195.68 * .3048; % wingpasn | m
S = 5500 * .3048^2; % wing planform area | m^2
Ix_b = 1.82*10^7*14.5939*(.3048)^2; % | kg m^2
Iz_b = 4.97*10^7*14.5939*(.3048)^2; % | kg m^2
Izx_b = 9.7*10^5*14.5939*(.3048)^2; % | kg m^2

% Calculated flight parameters
m = W/g; % 747 mass | kg
Ix = Ix_b*cos(xi)^2 + Iz_b*sin(xi)^2 + Izx_b*sin(2*xi); % stability frame | kg m^2
Iz = Iz_b*cos(xi)^2 + Ix_b*sin(xi)^2 - Izx_b*sin(2*xi); % stability frame | kg m^2
Izx = .5*(Iz_b - Ix_b)*sin(2*xi) - Izx_b*(sin(xi)^2 - cos(xi)^2); % stability frame | kg m^2
Ix_p = (Ix*Iz - Izx^2)/Iz; % | kg m^2
Iz_p = (Ix*Iz - Izx^2)/Ix; % | kg m^2
Izx_p = Izx/(Ix*Iz - Izx^2); % | kg m^2

% Nondimensional derivatives for Boeing 747
C_y = [-.8771, 0, 0];
C_l = [-.2797, -.3295, .304];
C_n = [.1946, -.04073, -.2737];

%% Problem 1
% Calculate stability derivatives
[y, l, n] = lat_derivs(C_y, C_l, C_n, rho, u0, b, S);
[Y, L, N] = lat_stability_transform(y, l, n, xi);

%% Problem 2
% Construct A matrix
A = [Y(1)/m, Y(2)/m, Y(3)/m - u0, g*cos(theta0); ...
    L(1)/Ix_p + Izx_p*N(1), L(2)/Ix_p + Izx_p*N(2), ...
    L(3)/Ix_p + Izx_p*N(3), 0; ...
    L(1)*Izx_p + N(1)/Iz_p, L(2)*Izx_p + N(2)/Iz_p, ...
    L(3)*Izx_p + N(3)/Iz_p, 0; ...
    0, 1, tan(theta0), 0];

%% Problem 3
% Find eigenvalues of A matrix
[evec, eval] = eig(A);
eval = diag(eval);

% eigenvalues
eval_d = eval(find(imag(eval), 1)); % Dutch roll mode
eval_s = eval(find(max(real(eval)) == real(eval), 1)); % spiral mode
eval_r = eval(find(min(real(eval)) == real(eval), 1)); % rolling mode

% time constants
tao_d = -1/real(eval_d);
tao_s = -1/eval_s;
tao_r = -1/eval_r;

% natural frequencies
omega_n_d = sqrt(real(eval_d)^2 + imag(eval_d)^2);

% damping ratios
zeta_d = -real(eval_d)/omega_n_d;

%% Problem 4
% Dutch roll approximation
y_v = A(1,1);
n_v = A(3,1);
n_r = A(3,3);

eval_d_app = roots([1, -(y_v + n_r), (y_v*n_r + u0*n_v)]);

%% Problems 5
% Lateral dynamics integration
tspan = [0 500];

% a) check equilibrium
y0_a = [0 0 0 0 0 0]';
[t_a, y_a] = ode45(@(t, y) lat_dynamics(t, y, A, u0, theta0), tspan, y0_a);

% b) perturbations from trim
% i. delta v = 10 m/s
y0_bi = [10 0 0 0 0 0]';
[t_bi, y_bi] = ode45(@(t, y) lat_dynamics(t, y, A, u0, theta0), tspan, y0_bi);

% ii. delta p = .15 rad/s
y0_bii = [0 .15 0 0 0 0]';
[t_bii, y_bii] = ode45(@(t, y) lat_dynamics(t, y, A, u0, theta0), tspan, y0_bii);

% iii. 
y0_biii = [-1.8563, -.4185, .0311, .6148, 0, 0]';
[t_biii, y_biii] = ode45(@(t, y) lat_dynamics(t, y, A, u0, theta0), tspan, y0_biii);

% iv.
y0_biv = [2.9477, -.0063, .0758, 1.2431, 0, 0]';
[t_biv, y_biv] = ode45(@(t, y) lat_dynamics(t, y, A, u0, theta0), tspan, y0_biv);

%% Printing
% 1
deriv_table = array2table([Y' L' N'], 'VariableNames', {'Y', 'L', 'N'});
fprintf('Problem 1\n');
disp(deriv_table);

% 2
fprintf('Problem 2\nA:\n');
disp(A);

% 4
fprintf('Problem 4\n');
fprintf('Dutch roll eigenvalue: %.4f + %.4fi\n', real(eval_d), imag(eval_d));
fprintf('Dutch roll approximation: %.4f + %.4fi\n', real(eval_d_app(1)), imag(eval_d_app(1)));

%% Plotting
% 3a
figure();
subplot(2,1,1);
plot(t_a, y_a(:,1));
xlabel('Time (s)');
ylabel('Y Velocity (m/s)');
title('Side Velocity vs Time in Trim');
subplot(2,1,2);
plot(t_a, y_a(:,2));
xlabel('Time (s)');
ylabel('Roll Rate (rad/s)');
title('Roll Rate vs Time in Trim');

figure();
subplot(2,1,1);
plot(t_a, y_a(:,3));
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
title('Yaw Rate vs Time in Trim');
subplot(2,1,2);
plot(t_a, y_a(:,4));
xlabel('Time (s)');
ylabel('Bank Angle (rad)');
title('Bank Angle vs Time in Trim');

figure();
subplot(2,1,1);
plot(t_a, y_a(:,5));
xlabel('Time (s)');
ylabel('Azimuth (rad)');
title('Azimuth vs Time in Trim');
subplot(2,1,2);
plot(t_a, y_a(:,6));
xlabel('Time (s)');
ylabel('Y Displacement (m)');
title('Displacement vs Time in Trim');

% 3b
% i.
figure();
subplot(2,1,1);
plot(t_bi, y_bi(:,1));
xlabel('Time (s)');
ylabel('Y Velocity (m/s)');
title('Side Velocity vs Time with Velocity Perturbation');
subplot(2,1,2);
plot(t_bi, y_bi(:,2));
xlabel('Time (s)');
ylabel('Roll Rate (rad/s)');
title('Roll Rate vs Time with Velocity Perturbation');

figure();
subplot(2,1,1);
plot(t_bi, y_bi(:,3));
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
title('Yaw Rate vs Time with Velocity Perturbation');
subplot(2,1,2);
plot(t_bi, y_bi(:,4));
xlabel('Time (s)');
ylabel('Bank Angle (rad)');
title('Bank Angle vs Time with Velocity Perturbation');

figure();
subplot(2,1,1);
plot(t_bi, y_bi(:,5));
xlabel('Time (s)');
ylabel('Azimuth (rad)');
title('Azimuth vs Time with Velocity Perturbation');
subplot(2,1,2);
plot(t_bi, y_bi(:,6));
xlabel('Time (s)');
ylabel('Y Displacement (m)');
title('Displacement vs Time with Velocity Perturbation');

% ii.
figure();
subplot(2,1,1);
plot(t_bii, y_bii(:,1));
xlabel('Time (s)');
ylabel('Y Velocity (m/s)');
title('Side Velocity vs Time with Roll Rate Perturbation');
subplot(2,1,2);
plot(t_bii, y_bii(:,2));
xlabel('Time (s)');
ylabel('Roll Rate (rad/s)');
title('Roll Rate vs Time with Roll Rate Perturbation');

figure();
subplot(2,1,1);
plot(t_bii, y_bii(:,3));
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
title('Yaw Rate vs Time with Roll Rate Perturbation');
subplot(2,1,2);
plot(t_bii, y_bii(:,4));
xlabel('Time (s)');
ylabel('Bank Angle (rad)');
title('Bank Angle vs Time with Roll Rate Perturbation');

figure();
subplot(2,1,1);
plot(t_bii, y_bii(:,5));
xlabel('Time (s)');
ylabel('Azimuth (rad)');
title('Azimuth vs Time with Roll Rate Perturbation');
subplot(2,1,2);
plot(t_bii, y_bii(:,6));
xlabel('Time (s)');
ylabel('Y Displacement (m)');
title('Displacement vs Time with Roll Rate Perturbation');

% iii.
figure();
subplot(2,1,1);
plot(t_biii, y_biii(:,1));
xlabel('Time (s)');
ylabel('Y Velocity (m/s)');
title('Side Velocity vs Time with Perturbation Set 1');
subplot(2,1,2);
plot(t_biii, y_biii(:,2));
xlabel('Time (s)');
ylabel('Roll Rate (rad/s)');
title('Roll Rate vs Time with Perturbation Set 1');

figure();
subplot(2,1,1);
plot(t_biii, y_biii(:,3));
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
title('Yaw Rate vs Time with Perturbation Set 1');
subplot(2,1,2);
plot(t_biii, y_biii(:,4));
xlabel('Time (s)');
ylabel('Bank Angle (rad)');
title('Bank Angle vs Time with Perturbation Set 1');

figure();
subplot(2,1,1);
plot(t_biii, y_biii(:,5));
xlabel('Time (s)');
ylabel('Azimuth (rad)');
title('Azimuth vs Time with Perturbation Set 1');
subplot(2,1,2);
plot(t_biii, y_biii(:,6));
xlabel('Time (s)');
ylabel('Y Displacement (m)');
title('Displacement vs Time with Perturbation Set 1');

% iv.
figure();
subplot(2,1,1);
plot(t_biv, y_biv(:,1));
xlabel('Time (s)');
ylabel('Y Velocity (m/s)');
title('Side Velocity vs Time with Perturbation Set 2');
subplot(2,1,2);
plot(t_biv, y_biv(:,2));
xlabel('Time (s)');
ylabel('Roll Rate (rad/s)');
title('Roll Rate vs Time with Perturbation Set 2');

figure();
subplot(2,1,1);
plot(t_biv, y_biv(:,3));
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
title('Yaw Rate vs Time with Perturbation Set 2');
subplot(2,1,2);
plot(t_biv, y_biv(:,4));
xlabel('Time (s)');
ylabel('Bank Angle (rad)');
title('Bank Angle vs Time with Perturbation Set 2');

figure();
subplot(2,1,1);
plot(t_biv, y_biv(:,5));
xlabel('Time (s)');
ylabel('Azimuth (rad)');
title('Azimuth vs Time with Perturbation Set 2');
subplot(2,1,2);
plot(t_biv, y_biv(:,6));
xlabel('Time (s)');
ylabel('Y Displacement (m)');
title('Displacement vs Time with Perturbation Set 2');