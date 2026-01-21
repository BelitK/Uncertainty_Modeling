% Make sure params exist
params_battery;   % creates struct p in workspace

% Operating point
xop = [0.4; 0; 0];
uop = 0;

% Finite difference step sizes
eps_x = 1e-6;     % for states
eps_u = 1e-6;     % for input

n = numel(xop);
m = 1;            % single input i
ny = 1;           % single output vT

% Evaluate at operating point
f0 = battery_f(xop, uop, p);
y0 = battery_h(xop, uop, p);

% Allocate
A = zeros(n,n);
B = zeros(n,m);
C = zeros(ny,n);
D = zeros(ny,m);

% A = df/dx, C = dh/dx
for k = 1:n
    dx = zeros(n,1);
    dx(k) = eps_x;

    fp = battery_f(xop + dx, uop, p);
    fm = battery_f(xop - dx, uop, p);
    A(:,k) = (fp - fm) / (2*eps_x);

    yp = battery_h(xop + dx, uop, p);
    ym = battery_h(xop - dx, uop, p);
    C(:,k) = (yp - ym) / (2*eps_x);
end

% B = df/du, D = dh/du
up = uop + eps_u;
um = uop - eps_u;

Bp = battery_f(xop, up, p);
Bm = battery_f(xop, um, p);
B(:,1) = (Bp - Bm) / (2*eps_u);

Dp = battery_h(xop, up, p);
Dm = battery_h(xop, um, p);
D(:,1) = (Dp - Dm) / (2*eps_u);

% Show results
disp("A="); disp(A);
disp("B="); disp(B);
disp("C="); disp(C);
disp("D="); disp(D);

% Save for Simulink
lin = struct();
lin.xop = xop;
lin.uop = uop;
lin.yop = y0;
lin.A = A; lin.B = B; lin.C = C; lin.D = D;

save("lin_battery.mat", "lin");
