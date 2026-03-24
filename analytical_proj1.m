tic

% Given Parameters
gam = 1.4;
R = 287; % J/kg-K
k = 2.63*10^-2; % W/m-K
mu = 1.79*10^-5; % N-s/m^2
Pr = 0.7;
t = 6.1 * 10^-3; % s 

p_1 = 10^4; % Pa
p_4 = 10^5; % Pa

p_r = 3.031; % p_2 / p_1

rho_1 = 0.125; % kg/m^3
rho_4 = 1; % kg/m^3

u_1 = 0; % m/s
u_4 = 0; % m/s

L_x = 10; % m
x_mid = 5; %m m
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
x_shock = x_mid + W * t;
x_left = x_mid - a_4 * t;
x_right = x_mid + (u_3 - a_3) * t;
x_contact = x_mid + u_2 * t;

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
        u(i) = 2 / (gam + 1) * (a_4 + (x(i) - x_mid) / t);
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

% Plotting
subplot(2,2,1); 
plot(x, P); 
xlabel('x-coordinate (m)'); 
ylabel('Pressure (Pa)'); grid on;

subplot(2,2,2); 
plot(x, rho); 
xlabel('x-coordinate (m)'); 
ylabel('Density (kg/m^3)'); 
grid on;

subplot(2,2,3); 
plot(x, u); 
xlabel('x-coordinate (m)'); 
ylabel('Velocity (m/s)'); 
grid on;

subplot(2,2,4); 
plot(x, T); 
xlabel('x-coordinate (m)'); 
ylabel('Temperature (K)'); 
grid on;

sgtitle(sprintf('Shock Tube Flow \n(Solution at t = %.1f ms)', t * 1000))

toc