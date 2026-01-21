function [xhat, yhat] = EKF(i, y)
% Discrete EKF for x = [sigma; vTS; vTL]
% Measurement: terminal voltage y = h(x,i) + v
%
% Requires in base workspace:
%   theta  (parameter vector)
% And functions on path:
%   battery_f(x,i,theta)
%   battery_h(x,i,theta)

% ----- tuning (simple, stable defaults) -----
Ts = 0.1;

% Process noise covariance (small)
Q = diag([1e-10, 1e-8, 1e-8]);

% Measurement noise variance: (10 mV)^2
R = 1e-4;

% ----- persistent filter state -----
persistent x P
if isempty(x)
    x = [0.4; 0; 0];        % initial estimate
    P = diag([1e-4, 1e-3, 1e-3]);  % initial covariance
end

% ----- prediction step -----
% Nonlinear dynamics (Euler discretization)
f = battery_f(x, i, theta);         % 3x1
x_pred = x + Ts * f;

% Clamp SOC for safety
x_pred(1) = min(max(x_pred(1), 0), 1);

% Jacobian A = df/dx via finite differences around x
A = jacobian_f_fd(x, i);

% Discrete-time linearization: F = I + Ts*A
F = eye(3) + Ts * A;

P_pred = F*P*F' + Q;

% ----- measurement prediction -----
yhat = battery_h(x_pred, i, theta);

% Jacobian C = dh/dx via finite differences around x_pred
C = jacobian_h_fd(x_pred, i);   % 1x3 row

% ----- update step -----
S = C*P_pred*C' + R;            % scalar
K = (P_pred*C') / S;            % 3x1

innov = y - yhat;
x = x_pred + K * innov;

% Clamp SOC after update
x(1) = min(max(x(1), 0), 1);

P = (eye(3) - K*C) * P_pred;

xhat = x;

end

% ===== helper: finite-difference Jacobians =====
function A = jacobian_f_fd(x0, i0)
epsx = 1e-6;
A = zeros(3,3);
f0 = battery_f(x0, i0, theta);
for k = 1:3
    x1 = x0;
    x1(k) = x1(k) + epsx;
    f1 = battery_f(x1, i0, theta);
    A(:,k) = (f1 - f0)/epsx;
end
end

function C = jacobian_h_fd(x0, i0)
epsx = 1e-6;
C = zeros(1,3);
y0 = battery_h(x0, i0, theta);
for k = 1:3
    x1 = x0;
    x1(k) = x1(k) + epsx;
    y1 = battery_h(x1, i0, theta);
    C(1,k) = (y1 - y0)/epsx;
end
end
