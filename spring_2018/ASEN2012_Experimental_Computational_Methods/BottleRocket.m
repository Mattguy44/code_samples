% ASEN 2012
% Project 2: Bottle Rocket Design
%
% ID: 213
% Created: 11/21/2017
% Last Modified: 12/7/2017

% Assume: Rocket always pointed in direction of velocity
%         Standard Atmosphere
%         No wind

% Purpose: Model all three stages of water bottle rocket flight

% Inputs: initial air pressure, initial volume of water, drag coefficient,
%           launch angle

% Outputs: results from numeerical integration

function [t,y] = BottleRocket(P_a0, V_w0, C_d, theta0)
% Constants
g = 9.81; % acceleration due to gravity, m/s^2
C_dis = 0.8; % discharge coefficient
rho_atm = 0.961; % atmospheric density, kg/m^3
V_b = 0.002; % volume of bottle, m^3
P_atm = 12.1 * 6894.76; % atmospheric pressure, Pa
gamma = 1.4; % specific heat ratio
rho_w = 1000; % density of water, kg/m^3
d_t = 0.021; % diameter of throat, m
d_b = 0.105; % diameter of bottle, m
R = 287; % gas constant of air, J/kg/K or m^2/s^2/K
m_b = 0.15; % mass of bottle w/ cone and fins, kg
T_a0 = 300; % air temperature in bottle, K
L_s = 0.5; % length of test stand, m
L_s_x = L_s*cos(theta0); % length of test stand in x-direction, m

% Derived constants
A_b = pi*d_b^2/4; % cross-sectional area of bottle, m^2
A_t = pi*d_t^2/4; % area of the throat, m^2
V_a0 = V_b - V_w0; % initial volume of air, m^3
m_a0 = P_a0*V_a0/R/T_a0; % initial mass of air, kg
m_w0 = rho_w*V_w0; % initial mass of water, kg
P_end = P_a0*(V_a0/V_b)^gamma; % pressure of air when water is gone, Pa

% Other initial conditions
x0 = 0; % initial horizontal dispalcement, m
y0 = 0.25; % initial vertical displacement, m
v_x0 = 0; % initial horizontal velocity, m/s
v_y0 = 0; % initial vertical velocity, m/s

% ODEs
% Phase 1
tspan1 = [0 5];
init1 = [x0 y0 v_x0 v_y0 m_w0 m_a0 V_a0];
opts1 = odeset('Events', @waterOut);
[t1,y1] = ode45(@flight1, tspan1, init1, opts1);
% Phase 2
tspan2 = [t1(end) 5];
init2 = [y1(end,1:4) 0 y1(end,6:7)];
opts2 = odeset('Events', @airOut, 'RelTol', 1e-3, 'AbsTol', 1e-3);
[t2,y2] = ode45(@flight2, tspan2, init2, opts2);
% Phase 3
tspan3 = [t2(end) 7];
init3 = y2(end,:);
[t3,y3] = ode45(@flight3, tspan3, init3);
% Full flight
t = [t1;t2;t3];
y = [y1;y2;y3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value, isterminal, direction] = waterOut(t, y)
    % Stops Phase 1 when air volume is equal to bottle volume
    
    value = y(7) - V_b;
    isterminal = 1;
    direction = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value, isterminal, direction] = airOut(t, y)
    % Stops Phase 2 when air pressure is equal to ambient pressure
    
    m_a = y(6);
    P = P_end*(m_a/m_a0)^gamma;

    value = P - P_atm;
    isterminal = 1;
    direction = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1
% water in bottle
function dt = flight1(t, y)
    % Calculates the derivatives of the input variables
    % y = [displacementX, displacementY, velocityX, velocityY, massWater, massAir, volume]
    d_x = y(1); % displacement in x-direction, m
    v_x = y(3); % velocity in x-direction, m/s
    v_y = y(4); % velocity in y-direction, m/s
    m_w = y(5); % mass of water, kg
    m_a = y(6); % mass of air, kg
    V_a = y(7); % volume of air, m^3
    
    % known rates of change
    dd_x = v_x; % rate of change of x-displacement, m/s
    dd_y = v_y; % rate of change of y-displacement, m/s
    dm_a = 0; % rate of change for air mass, kg/s
    
    % derived variables
    P = P_a0*(V_a0/V_a)^gamma; % air pressure, Pa or kg/m/s^2
    v = sqrt(v_x^2 + v_y^2); % rocket speed, m/s
    D = 0.5*rho_atm*v^2*C_d*A_b; % drag on rocket, N
    m = m_w + m_a + m_b; % mass of rocket, kg
    F = 2*C_dis*(P - P_atm)*A_t; % thrust, N
    v_e = sqrt(2*(P - P_atm)/rho_w); % water exit velocity, m/s

    % derived rates of change
    dm_w = -F/v_e; % mass flow rate of water, kg/s
    dV = C_dis*A_t*v_e; % air volume rate of change, m^3/s

    % Part 1.1
    % rocket on rail; assume this is while water is still in bottle.
    if (d_x < L_s_x)
        ang = theta0 - pi/64; % rocket angle with respect to ground, rad

        g_r = g*sin(ang); % gravitational acceleration along rail, m/s^2
        a = (F-D)/m - g_r; % rocket acceleration; should be - g_r
        dv_x = a*cos(ang); % rate of change of velocity in x-direction, m/s^2
        dv_y = a*sin(ang); % rate of change of velocity in y-direction, m/s^2
        
    % Part 1.2
    % rocket off rail; water still in bottle
    else
        ang = atan(v_y/v_x); % rocket angle with respect to ground, rad

        F_x = F*cos(ang); % thrust in x-direction, N
        D_x = D*cos(ang); % drag in x-direction, N
        dv_x = (F_x - D_x)/m; % rate of change of velocity in x-direction, m/s^2

        F_y = F*sin(ang); % thrust in y-direction, N
        D_y = D*sin(ang); % drag in y-direction, N
        dv_y = (F_y - D_y)/m - g; % rate of change of velocity in y-direction, m/s^2
    end
    
    % output derivatives
    dt = [dd_x, dd_y, dv_x, dv_y, dm_w, dm_a, dV]';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2
% only air in bottle
function dt = flight2(t, y)
    % Calculates the derivatives of the input variables
    % y = [displacementX, displacementY, velocityX, velocityY, massWater, massAir, volume]
    v_x = y(3); % velocity in x-direction, m/s
    v_y = y(4); % velocity in y-direction, m/s
    m_a = y(6); % mass of air, kg
    
    % known rates of change
    dd_x = v_x; % rate of change of x-displacement, m/s
    dd_y = v_y; % rate of change of y-displacement, m/s
    dm_w = 0; % rate of change for mass of water, kg/s
    
    % derived variables
    ang = atan(v_y/v_x); % rocket angle with respect to ground, rad
    rho_a = m_a/V_b; % density of air, kg/m^3
    P = P_end*(m_a/m_a0)^gamma; % inner air pressure, Pa
    P_cr = P*(2/(gamma+1))^(gamma/(gamma - 1)); % air pressure at Mach 1, Pa
    T = P/(rho_a*R); % air temperature, K
    v = sqrt(v_x^2 + v_y^2); % rocket speed, m/s
    D = 0.5*rho_atm*v^2*C_d*A_b; % drag on rocket, N
    m = m_a + m_b; % mass of rocket, kg
    
    % correction for over-stepping
    if(P < P_atm)
        dd_x = 0;
        dd_y = 0;
        dv_x = 0;
        dv_y = 0;
        dm_a = 0;
        dV = 0;
    else
        % air exiting bottle at Mach 1
        if (P_cr > P_atm)
            P_e = P_cr;
            T_e = T*(2/(gamma+1)); % temperature of discharging air, K
            v_e = sqrt(gamma*R*T_e); % air exit velocity, m/s
            
        % air exiting bottle below Mach 1
        else
            P_e = P_atm;
            M = sqrt(2/(gamma - 1)*((P/P_atm)^((gamma - 1)/gamma) - 1)); % air exit Mach number
            T_e = T/(1 + (gamma - 1)/2*M^2); % air exit temperature, K
            v_e = M*sqrt(gamma*R*T_e); % air exit velocity, m/s
        end

        % force calculation
        rho_e = P_e/(R*T_e); % air density at exit, kg/m^3
        dm_a = -C_dis*rho_e*A_t*v_e; % rate of change of air mass, kg/s
        F = -dm_a*v_e + (P_e - P_atm)*A_t; % thrust, N
        
        F_x = F*cos(ang); % thrust in x-direction, N
        F_y = F*sin(ang); % thrust in y-direction, N
        D_x = D*cos(ang); % drag in x-direction, N
        D_y = D*sin(ang); % drag in y-direction, N
        dv_x = (F_x - D_x)/m; % rate of change of velocity in x-direction, m/s^2
        dv_y = (F_y - D_y)/m - g; % rate of change of velocity in y-direction, m/s^2

        dV = 0; % air volume rate of change, m^3/s
    end

    % output derivatives
    dt = [dd_x, dd_y, dv_x, dv_y, dm_w, dm_a, dV]';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3
% freefall
function dt = flight3(t,y)
    % Calculates the derivatives of the input variables
    % y = [displacementX, displacementY, velocityX, velocityY, massWater, massAir, volume]
    d_y = y(2); % displacement in y-direction, m
    v_x = y(3); % velocity in x-direction, m/s
    v_y = y(4); % velocity in y-direction, m/s
    m_a = y(6); % mass of air, kg
    
    % known rates of change
    dm_a = 0; % air mass rate of change, kg/s
    dm_w = 0; % water mass rate of change, kg/s
    dV = 0; % air volume rate of chnage, m^3/s
    dd_x = v_x; % rate of change of x-displacement, m/s
    dd_y = v_y; % rate of change of y-displacement, m/s
    
    % derived variables
    ang = atan(v_y/v_x); % rocket angle with respect to ground, rad
    v = sqrt(v_x^2 + v_y^2); % rocket speed, m/s
    D = 0.5*rho_atm*v^2*C_d*A_b; % drag on rocket, N
    m = m_a + m_b; % mass of the rocket, kg

    D_x = D*cos(ang); % drag in x-direction, N
    D_y = D*sin(ang); % drag in y-direction, N

    dv_x = -D_x/m; % rate of change of velocity in x-direction, m/s^2
    dv_y = -D_y/m - g; % rate of change of velocity in y-direction, m/s^2

    % output derivatives
    dt = [dd_x, dd_y, dv_x, dv_y, dm_w, dm_a, dV]';
    
    % stops movement when the rocket hits the ground
    if (d_y < 0)
        dt = zeros(7,1);
    end
end

end % end of BottleRocket
