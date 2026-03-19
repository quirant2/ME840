% Project 1 - Shock Tube Case 1 Analytical Solution
clear all; clc;

% Initial Parameters for Case 1 [cite: 474-476]
P4 = 1e5; rho4 = 1.0; u4 = 0;
P1 = 1e4; rho1 = 0.125; u1 = 0;
gamma = 1.4; R = 287;
t = 0.0061; % Example time at 6.1 msec [cite: 650]
L = 10; x_mid = 5; 
Nx = 1000; x = linspace(0, L, Nx);

% Sound speeds [cite: 447, 610]
a1 = sqrt(gamma * P1 / rho1);
a4 = sqrt(gamma * P4 / rho4);

% Pressure Ratio PR = P2/P1 (given for Case 1) 
PR = 3.031;
P2 = PR * P1;
P3 = P2; % cite: 451, 607

% Calculate Region 2 Properties 
alpha = (gamma + 1) / (gamma - 1);
rho2 = rho1 * (1 + alpha * PR) / (alpha + PR);
u2 = u1 + a1 * (PR - 1) / sqrt(gamma * (gamma-1)/2 * (1 + alpha * PR));
u3 = u2; % cite: 451, 607

% Region 3 Density and Sound Speed [cite: 452, 615]
rho3 = rho4 * (P3 / P4)^(1/gamma);
a3 = sqrt(gamma * P3 / rho3);

% Shock Speed (W) [cite: 447, 625]
W = u1 + a1 * sqrt((gamma + 1)/(2 * gamma) * (PR - 1) + 1);

% Wave Locations [cite: 581, 626, 632]
x_shock = x_mid + W * t;
x_contact = x_mid + u2 * t;
x_fan_right = x_mid + (u3 - a3) * t;
x_fan_left = x_mid - a4 * t;

% Initialize vectors
P = zeros(1, Nx); rho = zeros(1, Nx); u = zeros(1, Nx);

for i = 1:Nx
    rel_x = x(i);
    if rel_x <= x_fan_left % Region 4 [cite: 578]
        P(i) = P4; rho(i) = rho4; u(i) = u4;
    elseif rel_x < x_fan_right % Expansion Fan [cite: 456, 627, 628]
        u(i) = 2/(gamma+1) * (a4 + (x(i)-x_mid)/t);
        rho(i) = rho4 * (1 - (gamma-1)/2 * (u(i)/a4))^(2/(gamma-1));
        P(i) = P4 * (rho(i)/rho4)^gamma;
    elseif rel_x <= x_contact % Region 3 [cite: 576]
        P(i) = P3; rho(i) = rho3; u(i) = u3;
    elseif rel_x <= x_shock % Region 2 [cite: 571, 572]
        P(i) = P2; rho(i) = rho2; u(i) = u2;
    else % Region 1 [cite: 570]
        P(i) = P1; rho(i) = rho1; u(i) = u1;
    end
end

% Plotting results [cite: 493, 498]
subplot(2,2,1); plot(x, P); title('Pressure (Pa)'); grid on;
subplot(2,2,2); plot(x, rho); title('Density (kg/m^3)'); grid on;
subplot(2,2,3); plot(x, u); title('Velocity (m/s)'); grid on;