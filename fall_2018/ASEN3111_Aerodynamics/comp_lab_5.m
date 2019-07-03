% ASEN 3111 - Aerodynamics
% Computational Lab 5
% CFD Simulation of Flow Over Aerodynamic Bodies
%
% Author: Matthew Ryan
% Created: Dec. 12, 2018
% Last Edit: Dec. 13, 2018

% Housekeeping
close all; clear; clc;

%% Read data
data = csvread('Fluent Data - Sheet1.csv',1,0);
aoa = data(:,1); % angle of attack | deg
aoa_rad = data(:,2); % angle of attack | rad
cl = data(:,5);
cd = data(:,6);

%% Additional data
[lift_slope, cl_aoa_0, aoa_L0, aoa_stall, cl_max, max_ind] = characteristics(aoa, cl);
cl_fit = aoa_rad*lift_slope + cl_aoa_0;

%% Plot data

% Sectional coefficient of lift vs angle of attack
figure();
hold on;
plot(aoa, cl, '.-');
plot(aoa, cl_fit);
xlabel(sprintf('Angle of Attack (%s)', char(176)));
ylabel('Sectional Coefficient of Lift');
title('c_l vs \alpha for NACA 0012 Airfoil');
legend('CFD Data', 'Best Fit Line', 'Location', 'best');
grid on;

% Drag polar
figure();
plot(cl,cd, '.-');
xlabel('Sectional Coefficient of Lift');
ylabel('Sectional Coefficient of Drag');
title('Drag Polar for NACA 0012 Airfoil (All Angles)');
grid on;

% drag polar limited to before stall
figure();
plot(cl(1:max_ind),cd(1:max_ind), '.-');
xlabel('Sectional Coefficient of Lift');
ylabel('Sectional Coefficient of Drag');
title('Drag Polar for NACA 0012 Airfoil (Before Stall)');
grid on;

%% Thin airfoil theory
lift_slope_TAT = 2*pi;
[~,ind_0] = min(abs(aoa));
aoa_L0_TAT = -cl(ind_0)/2/pi; % a_L0 = (a - c_l)/(2*pi)
% aoa_stall_TAT = N/A; % stall angle of attack | deg
% cl_max_TAT = 2*pi * aoa_stall_TAT*pi/180 + cl(ind_0); % max cl

%% Vortex panel method
% Airfoil Parameters
M = 0;
P = 0;
T = 12;
c = 1;
N = 300;

% Flow Parameters
v_inf = 1;
aoa_vortex_deg = -10:10:20;
aoa_vortex = aoa_vortex_deg*pi/180;
n_vortex = length(aoa_vortex);

[x,y] = NACA_Airfoils(M,P,T,c,N);

c_l_vortex = zeros(1, n_vortex);
for ii = 1:n_vortex
    c_l_vortex(ii) = Vortex_Panel(x, y, v_inf, aoa_vortex(ii));
end

% Find the lift slope and alpha_{L=0}
a_vortex = polyfit(aoa_vortex, c_l_vortex, 1);
lift_slope_vortex = a_vortex(1);
aoa_L0_vortex = -c_l_vortex(2)/lift_slope_vortex;

% Stall angle and cl_max
aoa_stall_vortex = 90;
cl_max_vortex = Vortex_Panel(x, y, v_inf, aoa_stall_vortex);

%% Experiment data
% Ladson 80 grit NACA 0012 at Re = 6000000
data_ladson_80 = dlmread('ladson_80grit.txt');
aoa_ladson_80 = data_ladson_80(:,1);
cl_ladson_80 = data_ladson_80(:,2);
cd_ladson_80 = data_ladson_80(:,3);

[lift_slope_ladson_80,~, aoa_L0_ladson_80, aoa_stall_ladson_80, cl_max_ladson_80, max_ind_ladson_80] = characteristics(aoa_ladson_80, cl_ladson_80);

% Ladson 120 grit NACA 0012 at Re = 6000000
data_ladson_120 = dlmread('ladson_120grit.txt');
aoa_ladson_120 = data_ladson_120(:,1);
cl_ladson_120 = data_ladson_120(:,2);
cd_ladson_120 = data_ladson_120(:,3);

[lift_slope_ladson_120,~, aoa_L0_ladson_120, aoa_stall_ladson_120, cl_max_ladson_120, max_ind_ladson_120] = characteristics(aoa_ladson_120, cl_ladson_120);

% Ladson 180 grit NACA 0012 at Re = 6000000
data_ladson_180 = dlmread('ladson_180grit.txt');
aoa_ladson_180 = data_ladson_180(:,1);
cl_ladson_180 = data_ladson_180(:,2);
cd_ladson_180 = data_ladson_180(:,3);

[lift_slope_ladson_180,~, aoa_L0_ladson_180, aoa_stall_ladson_180, cl_max_ladson_180, max_ind_ladson_180] = characteristics(aoa_ladson_180, cl_ladson_180);

%% Plot comparisons
figure();
hold on;
plot(aoa, cl, '.-'); % CFD
plot(aoa, lift_slope_TAT*aoa_rad, '.-'); % TAT
plot(aoa_vortex_deg, c_l_vortex, '.-'); % VPM
plot(aoa_ladson_80, cl_ladson_80, '.-'); % Ladson, 80 grit
plot(aoa_ladson_120, cl_ladson_120, '.-'); % Ladson, 120 grit
plot(aoa_ladson_180, cl_ladson_180, '.-'); % Ladson, 180 grit
xlabel(sprintf('Angle of Attack (%s)', char(176)));
ylabel('Sectional Coefficient of Lift');
title('c_l vs \alpha Comparison');
legend('CFD', 'Thin Airfoil', 'Vortex Panel', 'Ladson 80', 'Ladson 120', 'Ladson 180', 'Location', 'best');
grid on;

figure();
hold on;
plot(cl, cd, '.-'); % CFD
plot(cl_ladson_80, cd_ladson_80, '.-'); % Ladson, 80 grit
plot(cl_ladson_120, cd_ladson_120, '.-'); % Ladson, 120 grit
plot(cl_ladson_180, cd_ladson_180, '.-'); % Ladson, 180 grit
xlabel('Sectional Coefficient of Lift');
ylabel('Sectional Coefficient of Drag');
title('Drag Polar for NACA 0012 Airfoil (All Angles)');
legend('CFD', 'Ladson 80', 'Ladson 120', 'Ladson 180', 'Location', 'best');
grid on;

figure();
hold on;
plot(cl(1:max_ind), cd(1:max_ind), '.-'); % CFD
plot(cl_ladson_80(1:max_ind_ladson_80), cd_ladson_80(1:max_ind_ladson_80), '.-'); % Ladson, 80 grit
plot(cl_ladson_120(1:max_ind_ladson_120), cd_ladson_120(1:max_ind_ladson_120), '.-'); % Ladson, 120 grit
plot(cl_ladson_180(1:max_ind_ladson_180), cd_ladson_180(1:max_ind_ladson_180), '.-'); % Ladson, 180 grit
xlabel('Sectional Coefficient of Lift');
ylabel('Sectional Coefficient of Drag');
title('Drag Polar for NACA 0012 Airfoil (Before Stall)');
legend('CFD', 'Ladson 80', 'Ladson 120', 'Ladson 180', 'Location', 'best');
grid on;

%% Print results
fprintf('------------------------------------------\n');
fprintf('NACA 0012 CFD\n');
fprintf('------------------------------------------\n');
fprintf('Lift slope:\t\t\t%.4f cl/rad\n', lift_slope);
fprintf('Zero-lift angle of attack:\t%.4f%s\n', aoa_L0, char(176));
fprintf('Stall angle:\t\t\t%0.4f%s\n', aoa_stall, char(176));
fprintf('Maximum coefficient of lift:\t%.4f\n\n', cl_max);

fprintf('------------------------------------------\n');
fprintf('Thin Airfoil Theory\n');
fprintf('------------------------------------------\n');
fprintf('Lift slope:\t\t\t%.4f cl/rad\n', lift_slope_TAT);
fprintf('Zero-lift angle of attack:\t%.4f%s\n', aoa_L0_TAT, char(176));
fprintf('Stall angle:\t\t\t%s\n', 'N/A');
fprintf('Maximum coefficient of lift:\t%s\n\n', 'N/A');

fprintf('------------------------------------------\n');
fprintf('Vortex Panel Method\n');
fprintf('------------------------------------------\n');
fprintf('Lift slope:\t\t\t%.4f cl/rad\n', lift_slope_vortex);
fprintf('Zero-lift angle of attack:\t%.4f%s\n', aoa_L0_vortex, char(176));
fprintf('Stall angle:\t\t\t%0.4f%s\n', aoa_stall_vortex, char(176));
fprintf('Maximum coefficient of lift:\t%.4f\n\n', cl_max_vortex);

fprintf('------------------------------------------\n');
fprintf('Ladson 80 Grit Experiment\n');
fprintf('------------------------------------------\n');
fprintf('Lift slope:\t\t\t%.4f cl/rad\n', lift_slope_ladson_80);
fprintf('Zero-lift angle of attack:\t%.4f%s\n', aoa_L0_ladson_80, char(176));
fprintf('Stall angle:\t\t\t%0.4f%s\n', aoa_stall_ladson_80, char(176));
fprintf('Maximum coefficient of lift:\t%.4f\n\n', cl_max_ladson_80);

fprintf('------------------------------------------\n');
fprintf('Ladson 120 Grit Experiment\n');
fprintf('------------------------------------------\n');
fprintf('Lift slope:\t\t\t%.4f cl/rad\n', lift_slope_ladson_120);
fprintf('Zero-lift angle of attack:\t%.4f%s\n', aoa_L0_ladson_120, char(176));
fprintf('Stall angle:\t\t\t%0.4f%s\n', aoa_stall_ladson_120, char(176));
fprintf('Maximum coefficient of lift:\t%.4f\n\n', cl_max_ladson_120);

fprintf('------------------------------------------\n');
fprintf('Ladson 180 Grit Experiment\n');
fprintf('------------------------------------------\n');
fprintf('Lift slope:\t\t\t%.4f cl/rad\n', lift_slope_ladson_180);
fprintf('Zero-lift angle of attack:\t%.4f%s\n', aoa_L0_ladson_180, char(176));
fprintf('Stall angle:\t\t\t%0.4f%s\n', aoa_stall_ladson_180, char(176));
fprintf('Maximum coefficient of lift:\t%.4f\n\n', cl_max_ladson_180);

%% Functions
function [lift_slope, cl_aoa_0, aoa_L0, aoa_stall, cl_max, max_ind] = characteristics(aoa, cl)
    % Convert angle of attack to radians
    aoa_rad = aoa.*pi/180;
    
    % Maximum sectional coefficient of lift
    [cl_max,max_ind] = max(cl);

    % Stall angle
    aoa_stall = aoa(max_ind); % | deg

    % Lift slope
    % taken from -10 to 10 degrees, as suggested by McCroskey, W.J.
    [~,i_10] = min(abs(aoa - 10));
    [~,i_n10] = min(abs(aoa + 10));
    a_1 = polyfit(aoa_rad(i_n10:i_10), cl(i_n10:i_10), 1);
    cl_fit1 = aoa_rad.*a_1(1) + a_1(2); % fit from -10 to 10
    % remove outliers
    N = [];
    for n = 1:max_ind
        if(cl_fit1(n) - cl(n) < 0.01)
            N = [N, n];
        end
    end
    a = polyfit(aoa_rad(N), cl(N), 1); % [lift slope, cl at zero angle]
    lift_slope = a(1);
    cl_aoa_0 = a(2);

    % Zero lift angle of attack
    aoa_L0 = -a(2)/a(1) * 180/pi; % | deg
end