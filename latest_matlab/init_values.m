%% init_battery_params.m
% Parameter list theta (22x1). Replace placeholders with your real values.

iT=0.05;

CBat  = 3600;     % [Coulomb]

% R_S(soc) = a*exp(b*soc) + c
RS_a  = 0.01;  RS_b  = -2.0;  RS_c  = 0.005;

% R_TS(soc), C_TS(soc)
RTS_a = 0.02;  RTS_b = -2.0;  RTS_c = 0.01;
CTS_a = 2000;  CTS_b = -2.0;  CTS_c = 500;

% R_TL(soc), C_TL(soc)
RTL_a = 0.03;  RTL_b = -2.0;  RTL_c = 0.015;
CTL_a = 3000;  CTL_b = -2.0;  CTL_c = 800;

% OCV(soc) = v0*exp(v1*soc) + v2 + v3*soc + v4*soc^2 + v5*soc^3
v0 = 3.5;  v1 = -2.0;  v2 = 0.0;  v3 = 0.0;  v4 = 0.0;  v5 = 0.0;

theta = [ ...
    CBat; ...
    RS_a; RS_b; RS_c; ...
    RTS_a; RTS_b; RTS_c; ...
    CTS_a; CTS_b; CTS_c; ...
    RTL_a; RTL_b; RTL_c; ...
    CTL_a; CTL_b; CTL_c; ...
    v0; v1; v2; v3; v4; v5 ...
];

s0   = 0.4;
u0   = 0;
vTS0 = 0;
vTL0 = 0;

% ---- compute Jacobian once ----
make_jacobian_abcd(theta, s0, u0, vTS0, vTL0);