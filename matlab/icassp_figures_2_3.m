
clc, clear all, close all

% Parameters
K   = 5;
N   = 2 * K^2;
T   = 1 / 16;
tau = N*T;
P   = 15;
TTs = 64;
T_s = T / TTs;
n   = 1 : N + (P+1);
t_n = n*T;

% Construct E-spline kernel
m            = 0:P;
gamma        = (2*pi) / (N-P);
alpha_0      = -1j * gamma * P / 2;
alpha_m      = alpha_0 + 1j * m * gamma;
[phi, t_phi] = generate_e_spline(alpha_m, T_s, T, 'anticausal');
phi          = real(phi);

figure
set(gcf, 'Position', [100 100 280 210])
plot(-t_phi(end:-1:1), phi(end:-1:1), 'k')

% Construct Dirichlet kernel
t_h = (0 : T_s : (P+1)*T);
t_b = t_h / ((P+1)*T) * 2 * pi - pi;
b   = diric(t_b, (P+1));

figure
set(gcf, 'Position', [150 150 280 210])
plot(t_h, b, 'k')

% Construct E-spline kernel of different order
P            = 9;
m            = 0:P;
gamma        = (2*pi) / (N-P);
alpha_0      = -1j * gamma * P / 2;
alpha_m      = alpha_0 + 1j * m * gamma;
[phi, t_phi] = generate_e_spline(alpha_m, T_s, T, 'anticausal');
phi          = real(phi);

% Border effect in the left side, Dirac outside tau interval that affects
% samples y_n
t_1     = -2.5*T;
a_1     = 1;
n_1     = 1:N;
phi_mat = get_phi_tk_n_mat(phi, t_phi, t_1, n_1, T, T_s);
y_1     = phi_mat * a_1;

% Border effect in the right side, Dirac inside tau that affects samples
% after NT
t_2     = tau-3*T;
a_2     = 1;
phi_mat = get_phi_tk_n_mat(phi, t_phi, t_2, n, T, T_s);
y_2     = phi_mat * a_2;

figure
set(gcf, 'Position', [200 200 280 210])
stem(n_1*T, y_1, '.k')
hold on
stem(t_1, a_1, '^k', 'fill', 'markersize', 10)
stem([0 tau], [1.35 1.35], 'k', 'Marker', 'none')
plot([0 tau], [1.35 1.35], '+-k', 'LineWidth', 2)
text(tau/2, 1.5, '\tau', 'FontSize', 20);
axis([-(P+1)*T tau+(P+1)*T -.2 1.8])

figure
set(gcf, 'Position', [250 250 280 210])
stem(n*T, y_2, '.k')
hold on
stem(t_2, a_2, '^k', 'fill', 'markersize', 10)
stem([0 tau], [1.35 1.35], 'k', 'Marker', 'none')
plot([0 tau], [1.35 1.35], '+-k', 'LineWidth', 2)
text(tau/2, 1.5, '\tau', 'FontSize', 20);
axis([-(P+1)*T tau+(P+1)*T -.2 1.8])


