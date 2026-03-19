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
x_mid = 5;

T_1 = p_1 / (rho_1 * R);
T_4 = p_4 / (rho_4 * R);
a_1 = sqrt(gam * R * T_1);
a_4 = sqrt(gam * R * T_4);

% 
p_2 = p_r * p_1;
alpha = (gam + 1) / (gam - 1);
rho_2 = rho_1 * (1 + alpha * p_r) / (alpha + p_r);
u_2 = u_1 + (a_1 * (p_r - 1)/(1 + alpha * p_r) ^ (1/2) * 1 / sqrt(gam * (gam -1) / 2));
a_2 = a_1 * sqrt(p_r * (alpha + p_r) / (1 + alpha * p_r));
T_2 = p_2 / (rho_2 * R);

p_3 = p_2;
u_3 = u_2;

rho_3 = (p_3 / p_4 * rho_4 ^ gam) ^ (1 / gam);
a_3 = a_4 - (gam - 1) / 2 * u_3;

W = u_1 + (p_r - 1) * a_1 ^ 2 / (gam * (u_2 - u_1));
x_shock = x_mid + W * t;