clear

tic
% Setup
L_x = 10; 
N_x = 1001; 
dx = L_x/(N_x - 1);
x = 0:dx:L_x;
x_mid = 5;
CFL = 0.5;
C_x = 0.3;

% Fluid Properties
R = 287; Pr = 0.7; mu = 1.79 * 10^-5; gam = 1.4; k = 2.63 * 10^-2;
c_v = R/(gam - 1); c_p = gam * c_v;

% Initilization
rho = zeros(1, N_x); p = zeros(1, N_x); u = zeros(1, N_x); tau = zeros(1, N_x); q = zeros(1, N_x);

p_1 = 10^4; rho_1 = 0.125; u_1 = 0;
p_4 = 10^5; rho_4 = 1; u_4 = 0;
% T_1 = p_1 / (rho_1 * R); T_4 = p_4 / (rho_4 * R);

rho(x <= x_mid) = rho_4;
p(x <= x_mid) = p_4;
u(x <= x_mid) = u_4;

rho(x > x_mid) = rho_1;
p(x > x_mid) = p_1;
u(x > x_mid) = u_1;

% p_r = 3.031;

t = 0;
t_final = 0.0058;

T = p ./ (rho * R);
e = c_v * T;
U = [rho; rho .* u; rho .* e];
steps = 0;

while t < t_final
    rho = U(1,:);
    u = U(2,:) ./ rho;
    e = U(3,:) ./ rho;
    T = e ./ c_v;
    p = rho .* R .* T;
    a = sqrt(gam * R .* T);

    nu_prime = max(4/3 * mu./rho, mu./rho * gam/Pr);
    dt = min(CFL * ((abs(u)/dx) + (a/dx) + (2 * nu_prime/(dx^2))).^-1);
    if t + dt > t_final
        dt = t_final - t;
    end
    
    [F, J] = MC(U, p, u, T, mu, k, dx);
    S = visc(U, p, C_x);
    
    U_p = U;
    for i = 2:N_x-1
        U_p(:, i) = U(:, i) - dt * (xdb2(F, dx, i) - J(:, i)) + S(:, i);
    end
    
    rho_p = U_p(1,:);
    u_p = U_p(2,:) ./ rho_p;
    e_p = U_p(3,:) ./ rho_p;
    T_p = e_p ./ c_v;
    p_p = rho_p .* R .* T_p;
    
    [F_p, J_p] = MC(U_p, p_p, u_p, T_p, mu, k, dx);

    S_p = visc(U_p, p_p, C_x);
    
    U_c = U;
    for i = 2:N_x-1
        U_c(:, i) = (1/2) * (U(:, i) + U_p(:, i) - dt * (xdf2(F_p, dx, i) - J_p(:, i))) + S_p(:, i);
    end
    
    U = U_c;
    t = t + dt;
    steps = steps + 1;
end

rho_f = U(1,:);
u_f = U(2,:) ./ rho_f;
e_f = U(3,:) ./ rho_f;
T_f = e_f/ c_v;
p_f = rho_f .* R .* T_f;
x_f = x;

function db2_dx = xdb2(F, dx, i)
    db2_dx = (F(:,i) - F(:,i-1)) / dx;
end

function df2_dx = xdf2(F, dx, i)
    df2_dx = (F(:,i+1) - F(:,i)) / dx;
end

function dc2_dx = xdc2(F, dx, i)
    dc2_dx = (F(:,i+1) - F(:,i-1)) / (2*dx);
end

% function db4_dx = xdb4(F, dx, i)
%     db4_dx = (7*F(:,i) - 8*F(:,i-1) + F(:,i-2)) / (6*dx);
% end
% 
% function df4_dx = xdf4(F, dx, i)
%     df4_dx = (-7*F(:,i) + 8*F(:,i-1) - F(:,i-2)) / (6*dx);
% end

function [F, J] = MC(U, p, u, T, mu, k, dx)
    [~, N_x] = size(U);
    F = zeros(3, N_x); J = zeros(3, N_x);

    for i = 1:N_x
        if i > 1 && i < N_x
            tau = (4 * mu / 3) * (xdc2(u, dx, i));
            q = -k * xdc2(T, dx, i);
            du_dx = xdc2(u, dx, i);
        else
            tau = 0; q = 0; du_dx = 0;
        end

        F(1,i) = U(2, i);
        F(2,i) = U(2, i)*u(i) + p(i) - tau;
        F(3,i) = (U(3,i)*u(i)) + q;

        J(3, i) = (tau - p(i)) * du_dx;
    end
end

function S = visc(U, p, C_x)
    [~, N_x] = size(U);
    S = zeros(3, N_x);
    for i = 2:N_x-1
        S(:, i) = C_x * ((abs(p(i+1) - 2*p(i) + p(i-1))) / (p(i+1) + 2*p(i) + p(i-1))) * (U(:, i+1) - 2*U(:, i) + U(:, i-1));
    end
end

% Analytical Solution
% Given Parameters
gam = 1.4;
R = 287; % J/kg-K
k = 2.63*10^-2; % W/m-K
mu = 1.79*10^-5; % N-s/m^2
Pr = 0.7;
 

p_1 = 10^4; % Pa
p_4 = 10^5; % Pa

p_r = 3.031; % p_2 / p_1

rho_1 = 0.125; % kg/m^3
rho_4 = 1; % kg/m^3

u_1 = 0; % m/s
u_4 = 0; % m/s

L_x = 10; % m
x_mid = 5; % m
N = 1001;
x = linspace(0,L_x,N);

% Speed of sound for State 1 and 4
T_1 = p_1 / (rho_1 * R);
T_4 = p_4 / (rho_4 * R);
a_1 = sqrt(gam * R * T_1);
a_4 = sqrt(gam * R * T_4);

% State 2
p_2 = p_r * p_1;
alpha = (gam + 1) / (gam - 1);
rho_2 = rho_1 * (1 + alpha * p_r) / (alpha + p_r);
u_2 = u_1 + (a_1 * (p_r - 1)/(1 + alpha * p_r) ^ (1/2) * 1 / sqrt(gam * (gam -1) / 2));
a_2 = a_1 * sqrt(p_r * (alpha + p_r) / (1 + alpha * p_r));
T_2 = p_2 / (rho_2 * R);

% State 3
p_3 = p_2;
u_3 = u_2;
rho_3 = (p_3 / p_4 * rho_4 ^ gam) ^ (1 / gam);
T_3 = p_3 / (rho_3 * R);
a_3 = a_4 - (gam - 1) / 2 * u_3;

% Shock Speed
W = u_1 + (p_r - 1) * a_1 ^ 2 / (gam * (u_2 - u_1));

% Region Cutoffs
x_shock = x_mid + W * t_final;
x_left = x_mid - a_4 * t_final;
x_right = x_mid + (u_3 - a_3) * t_final;
x_contact = x_mid + u_2 * t_final;

% Initialize vectors
P = zeros(1,N); 
rho = zeros(1,N); 
u = zeros(1,N);
T = zeros(1,N);

for i = 1:N
    x_i = x(i);
    if x_i <= x_left % Region 4
        P(i) = p_4;
        rho(i) = rho_4;
        u(i) = u_4;
        T(i) = T_4;
    elseif x_i <= x_right % Expansion Fan Region
        u(i) = 2 / (gam + 1) * (a_4 + (x(i) - x_mid) / t_final);
        rho(i) = rho_4 * (1 - (gam - 1) / 2 * (u(i) / a_4)) ^ (2 / (gam - 1));
        T(i) = T_4 * (1 - (gam - 1) / 2 * (u(i) / a_4)) ^ 2;
        P(i) = rho(i) * R * T(i);
    elseif x_i <= x_contact % Region 3
        P(i) = p_3;
        rho(i) = rho_3;
        u(i) = u_3;
        T(i) = T_3;
    elseif x_i <= x_shock % Region 2
        P(i) = p_2;
        rho(i) = rho_2;
        u(i) = u_2;
        T(i) = T_2;
    else % Region 1
        P(i) = p_1;
        rho(i) = rho_1;
        u(i) = u_1;
        T(i) = T_1;
    end
end

colors = [0 0.4470 0.7410; ...      % Blue
          0.8500 0.3250 0.0980; ... % Red-Orange
          0.9290 0.6940 0.1250; ... % Yellow-Gold
          0.4940 0.1840 0.5560];    % Purple

figure('Name', 'MacCormack vs Exact Solution', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

% --- Density Plot ---
subplot(2,2,1); 
plot(x, rho, 'w-', 'LineWidth', 1.5); hold on;
plot(x_f, rho_f, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', colors(1,:), 'LineWidth', 0.7);
title('Density (\rho)', 'FontSize', 12); ylabel('kg/m^3'); grid on; box on;
legend('Exact','Numerical', 'Location', 'best');

% --- Velocity Plot ---
subplot(2,2,2); 
plot(x, u, 'w-', 'LineWidth', 1.5); hold on;
plot(x_f, u_f, 's', 'MarkerSize', 4, 'MarkerEdgeColor', colors(2,:), 'LineWidth', 0.7);
title('Velocity (u)', 'FontSize', 12); ylabel('m/s'); grid on; box on;

% --- Pressure Plot ---
subplot(2,2,3); 
plot(x, P, 'w-', 'LineWidth', 1.5); hold on;
plot(x_f, p_f, '^', 'MarkerSize', 4, 'MarkerEdgeColor', colors(3,:), 'LineWidth', 0.7);
title('Pressure (P)', 'FontSize', 12); ylabel('Pa'); xlabel('x (m)'); grid on; box on;

% --- Temperature Plot ---
subplot(2,2,4); 
plot(x, T, 'w-', 'LineWidth', 1.5); hold on;
plot(x_f, T_f, 'd', 'MarkerSize', 4, 'MarkerEdgeColor', colors(4,:), 'LineWidth', 0.7);
title('Temperature (T)', 'FontSize', 12); ylabel('K'); xlabel('x (m)'); grid on; box on;

% Add overall title with variable alignment
sgtitle(['Shock Tube Flow: Numerical vs. Exact Solutions (t = ', num2str(t_final*1000), ' ms)'], ...
    'FontSize', 14, 'FontWeight', 'bold');
toc