% ASEN 3200
% Oribtal Mechanics and Attitude Control
% Assignment 2
% Matthew Ryan
% Jan. 29, 2019

% Housekeeping
close all; clear all; clc;

% --------------------------------------------------------------------
%%%% [1] 9.7
% Given
m = [10 10 8 8 12 12]; % mass | kg
x = [1 -1 4 -2 3 -3]; % x position | m
y = [1 -1 -4 2 -3 3]; % y position | m
z = [1 -1 4 -2 -3 3]; % z position | m

% Center of mass
M = sum(m); % total mass | kg
x_cm = sum(m.*x)/M; % center of mass x | m
y_cm = sum(m.*y)/M; % center of mass y | m
z_cm = sum(m.*z)/M; % center of mass z | m

% Moment of inertia matrices
I_all = zeros(3,3,6);
for ii = 1:6
    x_rel = x(ii) - x_cm; % x position relative to cm | m
    y_rel = y(ii) - y_cm; % y position relative to cm | m
    z_rel = z(ii) - z_cm; % z position relative to cm | m
    m_ii = m(ii);
    I_all(1,1,ii) = m_ii*(y_rel^2 + z_rel^2);
    I_all(1,2,ii) = -m_ii*x_rel*y_rel;
    I_all(1,3,ii) = -m_ii*x_rel*z_rel;
    I_all(2,1,ii) = I_all(1,2,ii);
    I_all(2,2,ii) = m_ii*(x_rel^2 + z_rel^2);
    I_all(2,3,ii) = -m_ii*y_rel*z_rel;
    I_all(3,1,ii) = I_all(1,3,ii);
    I_all(3,2,ii) = I_all(2,3,ii);
    I_all(3,3,ii) = m_ii*(x_rel^2 + y_rel^2);
end
I_cm = sum(I_all,3); % moment of inertia matrix about cm | kg m^2

% Print
fprintf('[1] Curtis 9.7\n');
fprintf('I_cm:\n');
disp(I_cm);
fprintf('- - - - - - - - - - - - - - - - - - - - - -\n');


% --------------------------------------------------------------------
%%%% [2] 9.8
% Moment of inertia of the center of mass about the origin
I_cm_O = M*[y_cm^2 + z_cm^2,    -x_cm*y_cm,         -x_cm*z_cm; ...
            -y_cm*x_cm,         x_cm^2 + z_cm^2,    -y_cm*z_cm; ...
            -z_cm*x_cm,         -z_cm*y_cm,          x_cm^2 + y_cm^2];

% Moment of inertia about the origin
I_O = I_cm + I_cm_O;

% Moment of inertia about the axis through the origin
ax = [1 2 2];
ax_hat = ax/sqrt(ax*ax');
I_ax = ax_hat*I_O*ax_hat';

% Print
fprintf('[2] Curtis 9.8\n');
fprintf('I_ax:\n');
disp(I_ax);
fprintf('- - - - - - - - - - - - - - - - - - - - - -\n');


% --------------------------------------------------------------------
%%%% [3]
% Eigenvalues and eigenvectors
[I_evec, I_eval] = eig(I_cm);

% DCM
Q = [I_evec(:,1) I_evec(:,2) I_evec(:,3)]';

% Transformation
I_PA = Q*I_cm*Q'; % moment of inertia matrix about cm in principle axes

% Print
fprintf('[3]\n');
fprintf('Eigenvalues) %.2f, %.2f, %.2f\n', I_eval(1,1), I_eval(2,2), I_eval(3,3));
fprintf('Eigenvectors)\n');
fprintf('v_1:\n');
disp(I_evec(:,1));
fprintf('v_2:\n');
disp(I_evec(:,2));
fprintf('v_3:\n');
disp(I_evec(:,3));
fprintf('DCM to Principle Axis frame)\n');
fprintf('Q:\n');
disp(Q);
fprintf('Verification) Q*I_cm*Q'' = [v_1 v_2 v_3]\n');
fprintf('Q*I_cm*Q'':\n');
disp(I_PA);
fprintf('- - - - - - - - - - - - - - - - - - - - - -\n');


% --------------------------------------------------------------------
%%%% [4]
% Given
m_A = 3; % mass of cube | kg
m_B = 2; % mass of sphere | kg
l_A = 2; % length of cube side | m
r_B = 1; % radius of sphere | m
l_r = 2; % length of rod | m

% (a)
% Position of sections
y_A = -2; % location of cube cm | m
y_B = 2; % location of sphere cm | m

% Center of mass
x_cm4 = 0; % x position of cm | m
y_cm4 = (m_A*y_A + m_B*y_B)/(m_A + m_B); % y position of cm | m
z_cm4 = 0; % z position of cm | m

% Local moments of inertia
I_A = m_A*l_A^2/6; % moment of inertia of A about cm of A | kg m^2
I_B = .4*m_B*r_B^2; % moment of inertia of % about cm of B | kg m^2

% Moment of inertia about cm
I_cm4_x = I_A + m_A*(y_A - y_cm)^2 + I_B + m_B*(y_B - y_cm)^2;
I_cm4_y = I_A + I_B;
I_cm4_z = I_A + m_A*(y_A - y_cm)^2 + I_B + m_B*(y_B - y_cm)^2;
I_cm4 = diag([I_cm4_x, I_cm4_y, I_cm4_z]);

% (b)
% Given
w = 5; % spin rate | m/s
t = 3600; % time to stop | s

% Necessary torque
alpha = [0 0 -w/t]'; % angular acceleration | m/s^2
M_4 = I_cm4*alpha; 

% (c)
% Given
v_hat = [0 1/sqrt(2) 1/sqrt(2)];

% Moment of inertia about the origin
I_O4 = I_cm4 + diag([(m_A + m_B)*y_cm^2, 0, (m_A + m_B)*y_cm^2]);

% Moment of inertia about an axis through the origin
I_ax4 = v_hat*I_O4*v_hat';

% Print
fprintf('[4]\n');
fprintf('a) I_cm:\n');
disp(I_cm4);
fprintf('b) M:\n');
disp(M_4);
fprintf('c) I_ax:\n');
disp(I_ax4);
fprintf('- - - - - - - - - - - - - - - - - - - - - -\n');


% --------------------------------------------------------------------
%%%% [5]
% Given
m_5 = 1000; % mass of the spacecraft | kg
a_5 = [3, 1.5, 1]; % semi-major axes of effective ellipsoid | m

% (a) 
% Principle moments of inertia
I_5 = [.2*m_5*(a_5(2)^2 + a_5(3)^2), 0, 0; ...
       0, .2*m_5*(a_5(1)^2 + a_5(3)^2), 0; ...
       0, 0, .2*m_5*(a_5(1)^2 + a_5(2)^2)];

% (b)
% Given
w_5 = [2, -2, 2]'; % angular velocity | rad/s

% Angular Momentum
H_5 = I_5*w_5;
H_5_mag = sqrt(H_5'*H_5); % magnitude of angular momentum | kg m^2 / s

% Kinetic Energy
KE_5 = .5*w_5'*I_5*w_5; % | J

% (c)
% Counter-angular acceleration torques
% M = I d(w)/dt + w x (I w)
M_5 = cross(w_5, I_5*w_5);

% Print
fprintf('[5]\n');
fprintf('a) I_PA:\n');
disp(I_5);
fprintf('b)\nK = %.0f J\nH = %.0f kg m^2 / s\n\n', KE_5, H_5_mag);
fprintf('c) M:\n');
disp(M_5);
fprintf('- - - - - - - - - - - - - - - - - - - - - -\n');


% --------------------------------------------------------------------
%%%% [6]
% Given
m_6 = .15; % mass of a plate | kg
a_6 = .1; % side length of square plate | m

% Local moments of inertia
I_norm = m_6*a_6^2/6; % moment of inertia normal to plate | kg m^2
I_planar = m_6*a_6^2/12; % moment of inertia in plane of plate | kg m^2

% Moments of inertia about cm
I_norm_cm = I_norm;
I_planar_cm = I_planar + m_6*(.5*a_6)^2;

% Total moment of inertia matrix
I_6_PA = 2*I_norm_cm + 4*I_planar_cm; % principle axis moment of inertia
I_6 = diag([I_6_PA, I_6_PA, I_6_PA]);

% Print
fprintf('[6]\n');
fprintf('I_cube:\n');
disp(I_6);

% --------------------------------------------------------------------