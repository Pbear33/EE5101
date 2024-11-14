clc
clear
clear all
clf
close
close all

A = [-1.7, -0.25, 0; 23, -30, 20; 0, -940, -250];
B = [10, 0; -82, 0; 0, -670];
C = [0, 1, 0; 0, 0, 1];
D = 0;
x_0 = [1, 100, 200];


%design the semi-positive definite Q and the positive definite R
H = [10, 0, 0; 0, 10, 0];
Q = H'*H;
eig_q = eig(Q);
R = [1, 0; 0, 0.1];
eig_r = eig(R);

%A,B,H controllable and observable, guarantee the system stable 
Controllability = cat(2, B, A*B, A*A*B);
rank(Controllability)
Observability = cat(1, H, H*A, H*A*A);
rank(Observability)

% get the K, with ARE
gamma = cat(1, cat(2, A, -(B*inv(R))*B'), cat(2, -Q, -A'));

[v, d] = eig(gamma);
temp = cat(2, v(:, 1), v(:, 2), v(:, 3));
nu = cat(1, temp(1, :), temp(2, :), temp(3, :));
mu = cat(1, temp(4, :), temp(5, :), temp(6, :));

P = mu*inv(nu);
K = real(inv(R)*B'*P);




% observability
A1 = (A-B*K)';
W_o = cat(1, C, C*A1, C*A1*A1)
rank(W_o)

% choose the poles for oberver
desire_poles = [-5, -6, -7];
A_d = [0, 1, 0;0, 0, 1; -60, -47, -12];
L = place((A-B*K)', C', desire_poles)'

% perform pole placement
A1 = (A-B*K)';
B1 = C';
% first calculate the controllability matrix,W_c
W_c = cat(2, B1, A1*B1, A1*A1*B1)

% then group the independent vectors as Matrix_C, and its inverse
Matrix_C = cat(2, W_c(:, 1), W_c(:, 3), W_c(:, 2));
inv_Matrix_C = inv(Matrix_C);

% get the T, d1 = 2, d2 = 1
T = cat(1, inv_Matrix_C(2, :), inv_Matrix_C(2, :)*A1, inv_Matrix_C(3, :));

% get A_hat, B_hat
A_hat = T*A1/T;
B_hat = T*B1;

% cal K_hat and K
K_hat = pinv(B_hat)*(A_hat-A_d);
% L = (K_hat*T)';
% eig(A-B*K-L*C)






%simulate
dt = 0.0001;
T = 0.1;
n_steps = T/dt;

x = x_0';
x_hat = [0;0;0];
x_hist = zeros(3, n_steps);
x_hat_hist = zeros(3, n_steps);
u_hist = zeros(2, n_steps);

for k = 1:n_steps
    y = C * x;
    
    % LQR
    u = - K * x_hat; 
    
    %state update
    x_dot = A * x + B * u;
    x = x + x_dot * dt;
    
    % observer update
    x_hat_dot = A * x_hat + B * u + L * (y - C * x_hat);
    x_hat = x_hat + x_hat_dot * dt;
    
    % collect data
    x_hist(:, k) = x;
    x_hat_hist(:, k) = x_hat;
    u_hist(:, k) = u;
end


time = (0:n_steps-1) * dt;
figure;
subplot(3,1,1);
plot(time, x_hist(1,:), time, x_hat_hist(1,:), '--');
title('state x1');
legend('real', 'estimate');

subplot(3,1,2);
plot(time, x_hist(2,:), time, x_hat_hist(2,:), '--');
title('state x2');
legend('real', 'estimate');

subplot(3,1,3);
plot(time, x_hist(3,:), time, x_hat_hist(3,:), '--');
title('state x3');
legend('real', 'estimate');

