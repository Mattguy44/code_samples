% ASEN 2012
% Project 2: Bottle Rocket Design
%
% Purpose:  Recalculate the thrust at each time step
%           Produce a thrust profile with the inputted name
% Assumptions: Smooth transition from one thrust phase to the next
% Inputs: time, [x-displacement, y-displacement, x-velocity, y-velocity,
%               mass of water, mass of air, volume of air], initial
%               conditions, thrust profile plot title
% Outputs: None
%
% ID: 213
% Created: 12/5/2017
% Modified: 12/5/2017

function thrustProf(t, y, init, titl)
    P_a0 = init(1);
    V_w0 = init(2);

    % Constants
    C_dis = 0.8; % discharge coefficient
    V_b = 0.002; % volume of bottle, m^3
    P_atm = 12.1 * 6894.76; % atmospheric pressure, Pa
    gamma = 1.4; % specific heat ratio
    d_t = 0.021; % diameter of throat, m
    R = 287; % gas constant of air, J/kg/K or m^2/s^2/K
    T_a0 = 300; % air temperature in bottle, K

    % Useful constants
    A_t = pi*d_t^2/4; % area of the throat, m^2
    V_a0 = V_b - V_w0; % initial volume of air, m^3
    m_a0 = P_a0*V_a0/R/T_a0; % initial mass of air, kg
    P_end = P_a0*(V_a0/V_b)^gamma; % pressure of air when water is gone, Pa

    [N,~] = size(t);
    F = zeros(N,1);
    F = zeros(N,1);
    for i = 1:N
        % Calculates the thrust at each time
        % y = [displacementX, displacementY, velocityX, velocityY, massWater, massAir, volume]
        m_a = y(i,6); % mass of air, kg
        V_a = y(i,7); % volume of air, m^3

        P = P_a0*(V_a0/V_a)^gamma; % air pressure, Pa or kg/m/s^2
        F1(i) = 2*C_dis*(P - P_atm)*A_t; % thrust, N
        
        if (V_a < V_b)
            % Part 1
            % water in bottle
            
            F(i) = 2*C_dis*(P - P_atm)*A_t; % thrust, N
            endP1 = i;
            
        else
            % Part 2
            % only air in bottle
            
            rho_a = m_a/V_b; % density of air, kg/m^3
            P = P_end*(m_a/m_a0)^gamma; % inner air pressure, Pa
            P_cr = P*(2/(gamma+1))^(gamma/(gamma - 1)); % air pressure at Mach 1, Pa
            T = P/(rho_a*R);

            if(P > P_atm)
                if (P_cr > P_atm)
                    P_e = P_cr;
                    T_e = T*(2/(gamma+1)); % temperature of discharging air, K
                    v_e = sqrt(gamma*R*T_e); % air exit velocity, m/s
                else
                    P_e = P_atm;
                    M = sqrt(2/(gamma - 1)*((P/P_atm)^((gamma - 1)/gamma) - 1)); % air exit Mach number
                    T_e = T/(1 + (gamma - 1)/2*M^2); % air exit temperature, K
                    v_e = M*sqrt(gamma*R*T_e); % air exit velocity, m/s
                end

                rho_e = P_e/(R*T_e);
                dm_a = -C_dis*rho_e*A_t*v_e; % rate of change of air mass, kg/s
                F(i) = -dm_a*v_e + (P_e - P_atm)*A_t; % thrust, N
                endP2 = i;
                
            else
                % Part 3
                % freefall

                F(i) = 0;
                endP3 = i;
                
            end
        end
    end
    
    % Plot thrust profile
    figure
    hold on;
    plot(t(1:endP1+1), F1(1:endP1+1));
    plot([t(endP1+1); t(endP1+1:endP2)], [F1(endP1+1); F(endP1+1:endP2)]);
    plot(t(endP2:endP3), F(endP2:endP3));
    xlim([0 0.5]);
    title(titl);
    xlabel('Time (s)');
    ylabel('Thrust (N)');
    legend('Phase 1', 'Phase 2', 'Phase 3');
    
end