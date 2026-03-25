tic
% Setup
L_x = 10; 
N_x = 101; 
dx = L_x/(N_x - 1);
x = 0:dx:L_x;
x_mid = 5;
CFL = 0.5;

% Fluid Properties
R = 287; Pr = 0.7; mu = 1.79 * 10^-5; gam = 1.4; k = 2.63 * 10^-2;
c_v = R/(gam - 1); c_p = gam * c_v;

% Initilization
rho = zeros(1, N_x); p = zeros(1, N_x); u = zeros(1, N_x);

p_1 = 10^4; rho_1 = 0.125; u_1 = 0;
p_4 = 10^5; rho_4 = 1; u_4 = 0;

rho(x <= x_mid) = rho_4;
p(x <= x_mid) = p_4;
u(x <= x_mid) = u_4;

rho(x > x_mid) = rho_1;
p(x > x_mid) = p_1;
u(x > x_mid) = u_1;

p_r = 3.031;

T = p ./ (rho * R);
e = c_v * T;
U = [rho; rho .* u; rho .* e];

% MacCormak

toc