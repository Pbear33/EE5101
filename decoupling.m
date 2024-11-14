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
x_0 = [1; 100; 200];

% 
% temp = C(2, :)
% res = temp*B

% calculate the F and B, u = -Kx+Fr
B_star = cat(1, C(1, :)*B, C(2, :)*B);
C_star = cat(1, C(1, :)*A, C(2, :)*A)
F = inv(B_star)
K = inv(B_star) * C_star


% % the system is not stable, so try to add a pole placement 
% % selected poles
% A_d = [0, 1, 0;0, 0, 1; -0.08, -0.48, -1.4];
% 
% A1 = A - B*K
% B1 = B * F
% 
% % first calculate the controllability matrix,W_c
% W_c = cat(2, B1, A1*B1, A1*A1*B1);
% 
% % then group the independent vectors as Matrix_C, and its inverse
% Matrix_C = cat(2, W_c(:, 1), W_c(:, 3), W_c(:, 2));
% inv_Matrix_C = inv(Matrix_C);
% 
% % get the T, d1 = 2, d2 = 1
% T = cat(1, inv_Matrix_C(2, :), inv_Matrix_C(2, :)*A1, inv_Matrix_C(3, :));
% 
% % get A_hat, B_hat
% A_hat = T*A1/T;
% B_hat = T*B1;
% 
% % cal K_hat and K
% K_hat = pinv(B_hat)*(A_hat-A_d);
% K = K_hat*T




% simulate the system x_dot = (A-BK)x + BFr



runningtime = 20;
% Simulate the closed-loop response with initial condition x0
sys_cl = ss(A-B*K, B*F, C, 0);

t = linspace(0, runningtime, 1000);  % Simulation time
u = zeros(size(t));
v = ones(size(t));

[y, t_out, x] = lsim(sys_cl, [u;v], t);
y_tmax = max(y(:, 1));
y_t1max = max(y(:, 2));
y_min = min(y(:, 1))
y_t1min = min(y(:, 2))
y_t = y(1000, 1);
y_t1 = y(1000, 2);

% Plot response
figure;
% plot(t, initial_response);
plot(t_out, y);
xlabel('Time (s)');
ylabel('States');
% legend('C_a', 'T', 'T_j');
legend('T', 'T_j');
title('closed loop output')
grid on;

% Display the overshoot
overshoot = ((y_tmax - y_t)/y_t + (y_t1max- y_t1)/y_t1)/2
% Display the setting time
settime1 = 0;
settime12 = 0;
settime2 = 0;
settime22 = 0;
for n = 1:1000
    if y(n, 1) == y_min
        settime1 = n;
        break;
    end
end
for n = 1:1000
    if y(n, 1) >= (y_t*0.98) && n > settime1
        settime12 = n;
        break;
    end
end
for n = 1:1000
    if y(n, 2) == y_t1min
        settime2 = n;
        break;
    end
end
for n = 1:1000
    if y(n, 2) >= (y_t1*0.98) && n > settime2
        settime22 = n;
        break;
    end
end
settime = (settime22 - settime2 + settime12 - settime1)/2*(runningtime/1000)

% Display the eigenvalues of closed-loop system
closed_loop_poles = eig(A - B * K);
disp('Closed-loop poles:');
disp(closed_loop_poles);


% state response
sys_cl = ss(A - B * K, zeros(3,2), eye(3), zeros(3,2));

% Simulate response to initial condition
t = linspace(0, runningtime, 1000);  % Simulation time
u = -1*ones(size(t));

initial_response = initial(sys_cl, x_0, t);
% Plot response
figure;
plot(t, initial_response);
xlabel('Time (s)');
ylabel('States');
legend('C_a', 'T', 'T_j');
title('State Responses with State Feedback Control');
grid on;
