clc;
clear;
close all;

A = [-1.7, -0.25, 0; 23, -30, 20; 0, -940, -250];
B = [10, 0; -82, 0; 0, -670];
C = eye(3); 


desire_poles = [-1, -0.2+0.2j, -0.2-0.2j];
A_d = [0, 1, 0;0, 0, 1; -20, -29, -10];

% first calculate the controllability matrix,W_c
W_c = cat(2, B, A*B, A*A*B);

% then group the independent vectors as Matrix_C, and its inverse
Matrix_C = cat(2, W_c(:, 1), W_c(:, 3), W_c(:, 2));
inv_Matrix_C = inv(Matrix_C);

% get the T, d1 = 2, d2 = 1
T = cat(1, inv_Matrix_C(2, :), inv_Matrix_C(2, :)*A, inv_Matrix_C(3, :));

% get A_hat, B_hat
A_hat = T*A/T;
B_hat = T*B;

% cal K_hat and K
K_hat = pinv(B_hat)*(A_hat-A_d);
K1 = K_hat*T






% set point and original state
x_sp = [5; 250; 300];
x_0 = [1; 100; 200];

a = 7; b = 4; c = 2; d = 5; 
W = diag([a + b + 1, c + 4, d + 5]);
Q = W; 
R = [1, 0, 0; 0, 1, 0;0, 0, 1];     
% K = lqr(A-B*K1, B*K1, Q, R)


A1 = A-B*K1;
B1 = B*K1;

% get the K, with ARE
gamma = cat(1, cat(2, A1, -(B1*inv(R))*B1'), cat(2, -Q, -A1'));

[v, d] = eig(gamma);
temp = cat(2, v(:, 1), v(:, 3), v(:, 4));
nu = cat(1, temp(1, :), temp(2, :), temp(3, :));
mu = cat(1, temp(4, :), temp(5, :), temp(6, :));

P = mu*inv(nu);
K = real(inv(R)*B1'*P)





% SIMULATE
dt = 0.0001;
T = 0.1;
n_steps = T / dt;
x = x_0;
x_hist = zeros(3, n_steps);
u_hist = zeros(3, n_steps);


eig(A-B*K1-B*K1*K)




sys_cl = ss(A-B*K1-B*K1*K, B*K1*K, C, 0);

t = linspace(0, 1, 1000);  % Simulation time
u = cat(1, x_sp(1, 1)*ones(size(t)),x_sp(2, 1)*ones(size(t)),x_sp(3, 1)*ones(size(t)));

[y, t_out, x] = lsim(sys_cl, u, t);

% Plot response
figure;
% plot(t, initial_response);
plot(t_out, y);
xlabel('Time (s)');
ylabel('States');
legend('C_a', 'T', 'T_j');
title('closed loop output')
grid on;
