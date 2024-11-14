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

%design the semi-positive definite Q and the positive definite R
H = [1.7, 0, 0; 0, 1.7, 0];
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


runningtime = 0.2;
% Simulate the closed-loop response with initial condition x0
sys_cl = ss(A-B*K, B, C, 0);

t = linspace(0, runningtime, 1000);  % Simulation time
u = -1*ones(size(t));

[y, t_out, x] = lsim(sys_cl, [u;u], t);
y_tmax = max(y(:, 1));
y_t1max = max(y(:, 2));
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
    if y(n, 1) > (y_t*0.02)
        settime1 = n;
        break;
    end
end
for n = 1:1000
    if y(n, 1) > (y_t*0.98)
        settime12 = n;
        break;
    end
end
for n = 1:1000
    if y(n, 2) >= (y_t1*0.02)
        settime2 = n;
        break;
    end
end
for n = 1:1000
    if y(n, 2) >= (y_t1*0.98)
        settime22 = n - settime2;
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
