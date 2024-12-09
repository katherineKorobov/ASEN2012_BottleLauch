%{
Author: Katherine Korobov
Assignment: Project 2 Part 1
Purpose: Simulate a bottle rocket with water expulsion using physics
principles and ode45.
%}

clc;
clear;
close all;

%Set up Constants Structure
const = setConst();

%Initial Conditions
initialConditions.V0_air = const.VemptyB - const.V0_water; %[m^3] Initial Air Volume
initialConditions.m0_water = const.rhoWater * const.V0_water; %[kg] Initial Water Mass
initialConditions.m0_air = (const.p0 * initialConditions.V0_air) / (const.Rair * const.T0); %[kg] Ideal Gas Law
initialConditions.m0_rocket = const.mEmpty + initialConditions.m0_water + initialConditions.m0_air; %[kg] Initial Rocket Mass

%Define Initial Statevector
%statevector = [x, v_x, z, v_z, m_r, volume_air, m_air];
statevector_0 = [const.x0, const.v0, const.z0, const.v0, initialConditions.m0_rocket, initialConditions.V0_air, initialConditions.m0_air];

%Define Integration Time
final_time = 5;
tspan = [0, final_time];

[t,statevector] = ode45(@(t,statevector) bottleMotion(t,statevector, const, initialConditions), tspan, statevector_0);

%% Plot Results
%Initialize thrust array
thrust = zeros(size(statevector, 1), 1);

% Loop through each row of statevector
for i = 1:size(statevector, 1)
    % Calculate thrust using bottleMotion
    [~, thrust_temp] = bottleMotion(t(i), statevector(i, :), const, initialConditions);

    % Store thrust value (ensure thrust_temp is scalar)
    thrust(i) = norm(thrust_temp, 2);
end

load('project2verification.mat');

%Position Plot
figure();
hold on;
plot(statevector(:, 1), statevector(:, 3), 'c-');
plot(verification.distance, verification.height, 'k--');
legend('Simulated Code', 'Verification', 'Location', 'northwest');
title('Rocket Trajectory');
xlabel('Horizontal Position (m)');
ylabel('Vertical Position (m)');
grid on;
hold off;

%Thrust
figure();
hold on;
plot(t, thrust, 'c-');
plot(verification.time, verification.thrust, 'k--');
legend('Simulated Code', 'Verification', 'Location', 'best');
title('Thrust Over Time');
xlabel('Times (s)');
ylabel('Thrust (N)');
grid on;
xlim([0 0.2]);
hold off;

%% Saved Variables
 maxThrust = max(thrust); 
 maxHeight = max(statevector(:,3)); 
 maxDistance = max(statevector(:,1));

%% Phase function that calculate d_dt statevector to be put into ode45

function [d_statevector_dt, Fthrust] = bottleMotion(t, statevector, const, initialConditions)
    %Devectorize:
    x_pos = statevector(1);
    x_vel = statevector(2);
    z_pos = statevector(3);
    z_vel = statevector(4);
    m_rocket = statevector(5);
    V_air = statevector(6);
    m_air = statevector(7);

    %Compute Velocity Magnitude and Heading
    mag_vel = sqrt((x_vel)^2 + (z_vel)^2);

   if (x_pos < (const.standLength * cosd(const.Theta0)) || mag_vel == 0) %Either velocity is 0 or it is still on the stand
        h_hat = [cosd(const.Theta0), sind(const.Theta0)];
   else
        h_hat = [x_vel / mag_vel, z_vel / mag_vel];
   end

    exit_area = pi * (const.exitD / 2)^2; %Exit Area
    normal_area = pi * (const.bottleD / 2)^2; %Cross Sectional Area

    %Compute Air Pressure Using Isentropic Expansion
    p1_air = const.p0 * (initialConditions.V0_air / V_air)^const.gamma; %Equation 5

    %End Conditions of Air Expansion
    p_end = const.p0 * (initialConditions.V0_air / const.VemptyB)^const.gamma; %Equation 14
    p2_air = p_end * (m_air / initialConditions.m0_air)^const.gamma; %Equation 16 **

    if(V_air <= const.VemptyB) %Phase 1: Water is Expelled From Bottle
        %Calculate Exit Velocity
        exit_vel = sqrt((2*(p1_air - const.Pa)) / const.rhoWater); %Equation 8

        %Calculate Change in Mass of Water
        m_dotwater = const.cDis * const.rhoWater * exit_area * exit_vel; %Equation 6

        %Calculate Thrust
        Fthrust = h_hat * (m_dotwater * exit_vel); %Equation 9

        %Calculate Rate of Change of Air Volume
        Vdot_air = const.cDis * exit_area * exit_vel; %Equation 10
        
        mdot_rocket = -(m_dotwater); %Equation 11
        mdot_air = 0; %Mass of Air Remains Constant

    elseif((p2_air >= const.Pa)  && (V_air >= const.VemptyB)) %Phase 2: Air Expansion
        %Update Air Properties
        rho_air = m_air / const.VemptyB; %Equation 16
        T_air = p2_air / (rho_air * const.Rair); %Equation 16

        %Calcualte Critical Pressure
        p_crit = p2_air * (2/(const.gamma + 1))^(const.gamma/(const.gamma - 1)); %Equation 17
        
        if(p_crit > const.Pa) %If Flow is Choked
            exit_T = (2 / (const.gamma + 1)) * T_air; %Equation 18
            exit_p = p_crit; %Equation 18
            exit_rho = exit_p / (const.Rair * exit_T); %Equation 18
            exit_mach = 1;

        else %Where Critical Pressure < Ambient Pressure
            exit_mach = sqrt(2 * ((p2_air / const.Pa)^((const.gamma - 1) / const.gamma) - 1) / (const.gamma - 1)); % Equation 19
            exit_T = T_air / (1 + ((const.gamma - 1) / 2) * (exit_mach)^2); %Derived from Equation 19
            exit_p = const.Pa; %Equation 19
            exit_rho = const.Pa / (const.Rair * exit_T); %Equation 19
        end

        exit_vel = exit_mach * sqrt(const.gamma * const.Rair * exit_T); %Equation 18/19
        mdot_air = const.cDis * exit_rho * exit_area * exit_vel; %Equation 21
        Fthrust = h_hat * (mdot_air * exit_vel + (exit_p - const.Pa) * exit_area); 
        Vdot_air = 0;
        mdot_rocket = -(mdot_air);
        mdot_air = -(mdot_air); %Redfine for State Vector 

    else %Ballistic Phase
         Fthrust = 0;
         mdot_rocket = 0;
         Vdot_air = 0;
         mdot_air = 0;

         if(z_pos <= 0)
             d_statevector_dt = [0; 0; 0; 0; 0; 0; 0]; %Goes Below Ground, Nothing Changes
             return;
         end
    end

    %Compute Drag
    Fdrag = h_hat * (0.5 * const.rhoAir * (mag_vel)^2 * const.Cd * normal_area);
    
    %Compute Gravity
    Fgravity = [0, -(m_rocket * const.g)];

    %Compute Fnet
    Fnet = Fthrust - Fdrag + Fgravity;

    %Compute Acceleration
    a =  Fnet / m_rocket; 
    
    %Create Derivative of State Vector to be Put into ode45
    d_statevector_dt = [x_vel; a(1); z_vel; a(2); mdot_rocket; Vdot_air; mdot_air];
    
end

%% Function to Set Constants
function [const] = setConst()
    const.g = 9.81; %[m/s^2, acceleration due to gravity]
    const.cDis = 0.78; %[discharge coefficient]
    const.rhoAir = 0.961; %[kg/m^3, ambient air density]
    const.VemptyB = 0.002; %[m^3, volume of empty bottle]
    const.Pa = 12.1 * 6894.76; %[psia, atmospheric pressure, needs to be Pa]
    const.gamma = 1.4; %[ratio of specific heats for air]
    const.rhoWater = 1000; %[kg/m^3, water density]
    const.exitD = 2.1 / 100; %[m, diameter of the throat (exit) converted from cm]
    const.bottleD = 10.5 / 100; %[m, diameter of the bottle converted from cm]
    const.Rair = 287; %[J/(kg*K), specific gas constant of air]
    const.mEmpty = 0.15; %[kg, mass of empty 2-liter bottle with cone and fins]
    const.Cd = 0.425; %[drag coefficient]
    const.p0 = (48 + 12.1) * 6894.76; %[Pa, initial absolute pressure of air in bottle converted from psig]
    const.V0_water = 0.0005; %[m^3, initial volume of water inside bottle]
    const.T0 = 310; %[K, initial temperature of air]
    const.v0 = 0.0; %[m/s, initial velocity of rocket]
    const.Theta0 = 40; %[degree, initial angle of rocket]
    const.x0 = 0.0; %[m, initial horizontal distance]
    const.z0 = 0.25; %[m, initial vertical distance]
    const.standLength = 0.5; %[m, length of test stand]

end