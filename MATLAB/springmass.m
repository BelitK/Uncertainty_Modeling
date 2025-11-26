function exercise1_solver()
%% Params
p.m  = 1.0;
p.cF = 20.0;
p.bD = 2.0;
p.alpha = 5.0;
p.beta  = 1.0;

x0 = [0.1; 0.0];        % initial [z; zdot]
S0 = zeros(2,5);        % sensitivities wrt [m cF bD alpha beta]
tf = 10;                % sim time

%% Classical LQR on linearized model
[A,B] = linAB(p);
Qx = diag([10, 1]);
R  = 0.1;
Kx = lqr(A,B,Qx,R);

%% Sensitivity-aware LQR on augmented linear model
[A_aug,B_aug] = augAB(p);                     % 12x12 A_aug, 12x1 B_aug
Q_S = 1e2;                                    % weight for sensitivity states
Qaug = blkdiag(Qx, Q_S*eye(10));              % [x; vec(S)]
Kaug = lqr(A_aug, B_aug, Qaug, R);

%% Simulate nonlinear plant with both controllers
[t1,X1] = ode45(@(t,x) closed_loop_nl(t,x,@u_classic,p,Kx), [0 tf], x0);
[t2,X2] = ode45(@(t,xa) closed_loop_nl_aug(t,xa,@u_aug,p,Kaug), [0 tf], [x0; S0(:)]);

%% Plots
figure; grid on; hold on;
plot(t1, X1(:,1), 'LineWidth', 1.4);
plot(t2, X2(:,1), 'LineWidth', 1.4);
xlabel('t [s]'); ylabel('z(t)');
legend('Classical LQR','Sensitivity-aware LQR');

figure; grid on; hold on;
S_stack = X2(:,3:end);              % 10 sensitivity states after the 2 plant states
labels = {'m','cF','bD','alpha','beta'};
for i = 1:5
    idx = (i-1)*2 + (1:2);          % 2 rows per parameter -> [dz/dp_i; dzdot/dp_i]
    plot(t2, S_stack(:, idx(1)), 'LineWidth', 1.1); % show dz/dp_i
end
xlabel('t [s]'); ylabel('dz/dp_i');
legend(strcat('dz/d', labels));

disp('Done.');
end

%% Control laws
function u = u_classic(t, x, p, Kx)
% State feedback u = -Kx * x
u = -Kx * x;
end

function u = u_aug(t, xa, p, Kaug)
% Augmented feedback u = -Kaug * [x; vec(S)]
x = xa(1:2);
Svec = xa(3:end);
za = [x; Svec];
u = -Kaug * za;
end

%% Nonlinear closed loop with classical LQR
function dx = closed_loop_nl(~, x, ufun, p, Kx)
u = ufun([], x, p, Kx);
dx = f_nl(x, u, p);
end

%% Nonlinear closed loop with augmented controller
% Plant is nonlinear. We propagate linearized sensitivity states for control.
function dxa = closed_loop_nl_aug(~, xa, ufun, p, Kaug)
x = xa(1:2);
Svec = xa(3:end);
S = reshape(Svec, 2, 5);

u = ufun([], xa, p, Kaug);        % uses [x; vec(S)]
dx = f_nl(x, u, p);               % true nonlinear plant

% propagate sensitivities using the linearized form around the equilibrium
[A,B]     = linAB(p);
[dA, dB]  = dAB_dp(p);
dS = zeros(2,5);
for i = 1:5
    dS(:,i) = A*S(:,i) + dA{i}*x + dB{i}*u;  % dot S_i = A S_i + dA_i x + dB_i u
end

dxa = [dx; dS(:)];
end

%% Nonlinear dynamics: m zdd + bD zd + beta zd^3 + cF z + alpha z^3 = cF u
function dx = f_nl(x, u, p)
z = x(1); zd = x(2);
dx = [ zd;
      ( -p.bD*zd - p.beta*zd^3 - p.cF*z - p.alpha*z^3 + p.cF*u )/p.m ];
end

%% Linearization around x=0, u=0
% A = [0 1; -cF/m  -bD/m],  B = [0; cF/m]
function [A,B] = linAB(p)
A = [0 1;
    -p.cF/p.m  -p.bD/p.m];
B = [0; p.cF/p.m];
end

%% Partial derivatives of A and B wrt parameters [m cF bD alpha beta]
% alpha, beta do not appear in the equilibrium linearization
function [dA, dB] = dAB_dp(p)
% order: 1:m  2:cF  3:bD  4:alpha  5:beta
dA = cell(1,5); dB = cell(1,5);

dA{1} = [0 0;  p.cF/p.m^2   p.bD/p.m^2];   % dA/dm
dA{2} = [0 0; -1/p.m        0          ];   % dA/dcF
dA{3} = [0 0;  0           -1/p.m      ];   % dA/d bD
dA{4} = zeros(2);                           % dA/d alpha
dA{5} = zeros(2);                           % dA/d beta

dB{1} = [0; -p.cF/p.m^2];                   % dB/dm
dB{2} = [0;  1/p.m      ];                   % dB/dcF
dB{3} = [0;  0          ];                   % dB/dbD
dB{4} = [0;  0          ];                   % dB/d alpha
dB{5} = [0;  0          ];                   % dB/d beta
end

%% Build augmented LTI (x and sensitivities S stacked)
% State is za = [x; S1; S2; S3; S4; S5], with each S_i in R^2
% dot x   = A x + B u
% dot S_i = A S_i + dA_i x + dB_i u
function [A_aug, B_aug] = augAB(p)
[A,B] = linAB(p);
[dA, dB] = dAB_dp(p);

A_xx = A;
A_xS = zeros(2, 10);

A_Sx = zeros(10, 2);
for i = 1:5
    A_Sx( (i-1)*2+(1:2), : ) = dA{i};
end

A_SS = kron(eye(5), A);

A_aug = [A_xx, A_xS;
         A_Sx, A_SS];

B_S = zeros(10,1);
for i = 1:5
    B_S( (i-1)*2+(1:2) ) = dB{i};
end
B_aug = [B; B_S];
end
