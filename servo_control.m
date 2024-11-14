clc;
clear;
close all;

A = [-1.7, -0.25, 0; 23, -30, 20; 0, -940, -250];
B = [10, 0; -82, 0; 0, -670];
C = [0, 1, 0; 0, 0, 1];
x_0 = [1; 100; 200];  

% the augmented system
A_aug = [A, zeros(3, 2); -C, zeros(2, 2)];
B_aug = [B; zeros(2, 2)];
temp = cat(1, cat(2, A, B), cat(2, -C, zeros(2, 2)));
rank(temp);

% feedbak gain K and K_z
desired_poles_aug = [-2+1i, -2-1i, -1, -0.2+1i, -0.2-1i];  
A_d = [0, 1, 0, 0, 0; 0, 0, 1, 0, 0;-287.1474, 273.0073, -2.6539, 2.4782, -2.3915;0, 0, 0, 0, 1;286.7741, 315.6989, 115.1029, -2.4931, -2.7461];
eig(A_d)



% first calculate the controllability matrix,W_c
W_c = cat(2, B_aug, A_aug*B_aug, A_aug*A_aug*B_aug, A_aug*A_aug*A_aug*B_aug, A_aug*A_aug*A_aug*A_aug*B_aug);

% then group the independent vectors as Matrix_C, and its inverse
Matrix_C = cat(2, W_c(:, 1), W_c(:, 3), W_c(:, 5), W_c(:, 2), W_c(:, 4));
rank(Matrix_C)
inv_Matrix_C = inv(Matrix_C);

% get the T, d1 = 2, d2 = 1
T = cat(1, inv_Matrix_C(3, :), inv_Matrix_C(3, :)*A_aug, inv_Matrix_C(3, :)*A_aug*A_aug, inv_Matrix_C(5, :), inv_Matrix_C(5, :)*A_aug);

% get A_hat, B_hat
A_hat = T*A_aug/T;
B_hat = T*B_aug;

% cal K_hat and K
K_hat = pinv(B_hat)*(A_hat-A_d);
K = K_hat*T
K_x = K(:, 1:3);    
K_z = K(:, 4:5);    

% simulate
dt = 0.01;
T = 30;   
n_steps = T / dt;
r = [100; 150];  % setpoint


x = x_0;             
z = [0; 0];            
w = [0; 0];           
x_hist = zeros(3, n_steps);
y_hist = zeros(2, n_steps);
u_hist = zeros(2, n_steps);


for k = 1:n_steps
    if k * dt >= 10
        w = [-2; 5];  % disturbance
    end
    
  
    y = C * x;
    e = r - y;
    z = z + e * dt;
    u = -K_x * x - K_z * z;
    
    %update
    x_dot = A * x + B * u + [0; w];  
    x = x + x_dot * dt;
    
    x_hist(:, k) = x;
    y_hist(:, k) = y;
    u_hist(:, k) = u;
end


time = (0:n_steps-1) * dt;
% figure;
% subplot(3,1,1);
% plot(time, x_hist(1,:));
% title('state x1 (Ca)');
% xlabel('Time (s)');
% 
% subplot(3,1,2);
% plot(time, x_hist(2,:));
% title('state x2 (T)');
% xlabel('Time (s)');
% 
% subplot(3,1,3);
% plot(time, x_hist(3,:));
% title('state x3 (Tj)');
% xlabel('Time (s)');

figure;
plot(time, y_hist(1,:), time, y_hist(2,:));
title('output y');
legend('T', 'Tj');
xlabel('Time (s)');

figure;
plot(time, u_hist(1,:), time, u_hist(2,:));
title('input u');
legend('u1', 'u2');
xlabel('Time (s)');




