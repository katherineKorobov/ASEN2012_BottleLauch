%{
Author: Katherine Korobov, Alexander Lu
Assignment: Project 2 Part 2
Purpose: Simulate a bottle rocket with water expulsion using physics
principles and ode45. We will explore how different parameters affect
the trajectory and resulting thrust of the bottle rocket. Then, we will hit
a target 92m away with an range of +-1m.
%}

clc;
clear;
close all;

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
%% Phase function that calculate d_dt statevector to be put into ode45
function [d_statevector_dt, Fthrust, phaseNumber] = bottleMotion(t, statevector, const, initialConditions)
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

        phaseNumber = 1;

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

        phaseNumber = 2;

    else %Ballistic Phase
         Fthrust = 0;
         mdot_rocket = 0;
         Vdot_air = 0;
         mdot_air = 0;

         phaseNumber = 3;

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
%% Initial Conditions Function
function [initialConditions, statevector_0] = initializeVar(const)
    initialConditions.V0_air = const.VemptyB - const.V0_water; %[m^3] Initial Air Volume
    initialConditions.m0_water = const.rhoWater * const.V0_water; %[kg] Initial Water Mass
    initialConditions.m0_air = (const.p0 * initialConditions.V0_air) / (const.Rair * const.T0); %[kg] Ideal Gas Law
    initialConditions.m0_rocket = const.mEmpty + initialConditions.m0_water + initialConditions.m0_air; %[kg] Initial Rocket Mass

    %Define Initial Statevector: statevector = [x, v_x, z, v_z, m_r, volume_air, m_air];
    statevector_0 = [const.x0, const.v0, const.z0, const.v0, initialConditions.m0_rocket, initialConditions.V0_air, initialConditions.m0_air];
end
%% Thrust Calculation Function
function [thrust] = getThrust(t, statevector, const, initialConditions)
     thrust = zeros(length(t), 1);

     for j = 1:length(t)
        %Calculate thrust using bottleMotion
        [~, thrust_temp, ~] = bottleMotion(t(j), statevector(j, :), const, initialConditions);

        %Store thrust value (ensure thrust_temp is scalar)
        thrust(j) = norm(thrust_temp, 2);
     end
end
%% Initialize variables
const = setConst();

%Create Looping Vectors
allPressures = [35, 40, 45, 50, 55, 68]; %[psi] Units Converted Later
allV0Water = [0.0002, 0.0005, 0.0008]; %[m^3];
allCd = [0.3, 0.35, 0.40, 0.45, 0.5, 0.6];
allTheta0 = [10, 20, 40, 60, 80]; %[Degree]

%Define Integration Time
final_time = 10;
tspan = [0, final_time];

%% Explore Pressure
fig1 = figure('Name','Rocket Trajectory for Varying Initial Air Pressures');
hold on;
title('Rocket Trajectory for Varying Initial Air Pressures');
xlabel('Horizontal Position (m)');
ylabel('Vertical Position (m)');
ylim([0, 35]);
grid on;

fig2 = figure('Name','Thrust Over Time for Varying Initial Air Pressures');
hold on;
hold on;
title('Thrust Over Time for Varying Initial Air Pressures');
xlabel('Time (s)');
ylabel('Thrust (N)');
xlim([0 0.2]);
grid on;

for i = 1:length(allPressures)
    const.p0 = (allPressures(i) + 12.1) * 6894.76; %Convert to Pa
    [initialConditions, statevector_0] = initializeVar(const);

    %Simulate using ode45
    [t,statevector] = ode45(@(t,statevector) bottleMotion(t,statevector, const, initialConditions), tspan, statevector_0);

    figure(1);
    plot(statevector(:, 1), statevector(:, 3));
    
    thrust = getThrust(t, statevector, const, initialConditions);
    
    figure(2);
    plot(t, thrust);
end
figure(1);
legend('Pressure = 35psi', 'Pressure = 40psi', 'Pressure = 45psi', 'Pressure = 50psi', 'Pressure = 55psi', 'Pressure = 68psi', 'Location','northwest');
hold off;

figure(2);
legend('Pressure = 35psi', 'Pressure = 40psi', 'Pressure = 45psi', 'Pressure = 50psi', 'Pressure = 55psi', 'Pressure = 68psi', 'Location','northeast');
hold off;

%% Explore Initial Water Mass 
% Note: The density of water does not change, nor would it make sense to
% change it because it is a physical property of the material. Instead we
% will vary the values of the initial volume of the water inside the
% bottle.

fig3 = figure('Name','Rocket Trajectory for Varying Initial Water Masses');
hold on;
title('Rocket Trajectory for Varying Initial Water Masses');
xlabel('Horizontal Position (m)');
ylabel('Vertical Position (m)');
ylim([0, 35]);
grid on;

fig4 = figure('Name','Thrust Over Time for Varying Initial Water Masses');
hold on;
hold on;
title('Thrust Over Time for Varying Initial Water Masses');
xlabel('Time (s)');
ylabel('Thrust (N)');
xlim([0 0.2]);
grid on;

for i = 1:length(allV0Water)
    const.V0_water = allV0Water(i);
    [initialConditions, statevector_0] = initializeVar(const);

    %Simulate using ode45
    [t,statevector] = ode45(@(t,statevector) bottleMotion(t,statevector, const, initialConditions), tspan, statevector_0);

    figure(3);
    plot(statevector(:, 1), statevector(:, 3));
    
    thrust = getThrust(t, statevector, const, initialConditions);

    figure(4);
    plot(t, thrust);
end
figure(3);
legend('Volume = 0.0002', 'Volume = 0.0005', 'Volume = 0.0008', 'Location','northwest');
hold off;

figure(4);
legend('Volume = 0.0002', 'Volume = 0.0005', 'Volume = 0.0008','Location','northeast');
hold off;

%% Explore Coefficient of Drag
fig5 = figure('Name','Rocket Trajectory for Varying Coefficient of Drag');
hold on;
title('Rocket Trajectory for Varying Coefficient of Drag');
xlabel('Horizontal Position (m)');
ylabel('Vertical Position (m)');
ylim([0, 35]);
grid on;

fig6 = figure('Name','Thrust Over Time for Varying Coefficients of Drag');
hold on;
title('Thrust Over Time for Varying Coefficients of Drag');
xlabel('Time (s)');
ylabel('Thrust (N)');
xlim([0 0.2]);
grid on;

for i = 1:length(allCd)
    const.Cd = allCd(i);
    
    [initialConditions, statevector_0] = initializeVar(const);

    %Going to create a statevector
    [t,statevector] = ode45(@(t,statevector) bottleMotion(t,statevector, const, initialConditions), tspan, statevector_0);

    figure(5);
    plot(statevector(:, 1), statevector(:, 3));

     %Loop through each row of statevector
    thrust = getThrust(t, statevector, const, initialConditions);
    
    figure(6);
    plot(t, thrust);
end

figure(5);
legend('Cd = 0.3', 'Cd = 0.35', 'Cd = 0.4', 'Cd = 0.45', 'Cd = 0.5', 'Cd = 0.6', 'Location','northwest');
hold off;

figure(6);
legend('Cd = 0.3', 'Cd = 0.35', 'Cd = 0.4', 'Cd = 0.45', 'Cd = 0.5', 'Cd = 0.6', 'Location','northeast');
hold off;

%% Explore Launch Angle
fig7 = figure('Name','Rocket Trajectory for Varying Launch Angles');
hold on;
title('Rocket Trajectory for Varying Launch Angles');
xlabel('Horizontal Position (m)');
ylabel('Vertical Position (m)');
ylim([0, 50]);
grid on;

fig8 = figure('Name','Thrust Over Time for Varying Launch Angles');
hold on;
title('Thrust Over Time for Varying Launch Angles');
xlabel('Time (s)');
ylabel('Thrust (N)');
xlim([0 0.2]);
grid on;

for i = 1:length(allTheta0)
    const.Theta0 = allTheta0(i);

    [initialConditions, statevector_0] = initializeVar(const);

    %Going to create a statevector
    [t,statevector] = ode45(@(t,statevector) bottleMotion(t,statevector, const, initialConditions), tspan, statevector_0);

    figure(7);
    plot(statevector(:, 1), statevector(:, 3));

     %Loop through each row of statevector
    thrust = getThrust(t, statevector, const, initialConditions);
    
    figure(8);
    plot(t, thrust);
end

figure(7);
legend('Angle = 10', 'Angle = 20', 'Angle = 40', 'Angle = 60', 'Angle = 80', 'Location','northeast');
hold off;

figure(8);
legend('Angle = 10', 'Angle = 20', 'Angle = 40', 'Angle = 60', 'Angle = 80', 'Location','northeast');
hold off;

%% Interpolate for exact value to hit target
%Lets use a range of values and then interpolate to find an exact number
const = setConst();

pressureRange0 = linspace(50, 100, 100); %[psi]

required_distance = 92; %[m]

max_ranges = zeros(length(pressureRange0), 1);

for i = 1:length(pressureRange0)
    const.p0 = (pressureRange0(i) + 12.1) * 6894.76; %[Pa]

    [initialConditions, statevector_0] = initializeVar(const);

    % Simulate using ode45
    [~, statevector] = ode45(@(t, statevector) bottleMotion(t, statevector, const, initialConditions), tspan, statevector_0);

    % Record the maximum horizontal distance
    max_ranges(i) = max(statevector(:, 1));
end

% Interpolate to find the pressure that gives the target range
targetPressure = interp1(max_ranges, pressureRange0, required_distance, 'linear');


%% Hit Target
const = setConst();
const.p0 = (targetPressure + 12.1) * 6894.76;

[initialConditions, statevector_0] = initializeVar(const);

[t,statevector] = ode45(@(t,statevector) bottleMotion(t,statevector, const, initialConditions), tspan, statevector_0);

thrust = getThrust(t, statevector, const, initialConditions);

test_range = max(statevector(:, 1));

%get phases
phases = zeros(length(t), 1);
for i = 1:length(t)
    [~, ~, phases(i)] = bottleMotion(t(i), statevector(i, :), const, initialConditions);
end

phaseChange1  = find(phases == 2, 1);
phaseChange2 = find(phases == 3, 1);

phaseChange1Time = t(phaseChange1);
phaseChange2Time = t(phaseChange2);

figure('Name','Target Trajectory');
hold on;
plot(statevector(:, 1), statevector(:, 3));
xline(phaseChange1Time, 'r--'); % Phase 1 -> 2
xline(phaseChange2Time, 'r--'); % Phase 2 -> 3
title('Rocket Trajectory');
xlabel('Horizontal Position (m)');
ylabel('Vertical Position (m)');
grid on;
ylim([0 30]);
hold off;

%Thrust
figure('Name', 'Target Thrust');
hold on;
plot(t, thrust);
xline(phaseChange1Time, 'r--'); % Phase 1 -> 2
xline(phaseChange2Time, 'r--'); % Phase 2 -> 3
title('Thrust Over Time');
xlabel('Times (s)');
ylabel('Thrust (N)');
grid on;
xlim([0 0.2]);
hold off;