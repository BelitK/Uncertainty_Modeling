%% init_battery_params.m
% Parameter list theta (22x1). Replace placeholders with your real values.

iT=0.05;

CBat  = 3600;     % [Coulomb]

% R_S(soc) = a*exp(b*soc) + c
RS_a  = 0.25;  RS_b  = -20.0;  RS_c  = 0.07;

% R_TS(soc), C_TS(soc)
RTS_a = 1.0;  RTS_b = -30.0;  RTS_c = 0.015;
CTS_a = -900;  CTS_b = -2.0;  CTS_c = 1000;

% R_TL(soc), C_TL(soc)
RTL_a = 0.01;  RTL_b = -4.0;  RTL_c = 0.005;
CTL_a = 25000;  CTL_b = -2.0;  CTL_c = 2000;

% OCV(soc) = v0*exp(v1*soc) + v2 + v3*soc + v4*soc^2 + v5*soc^3
v0 = -1.0;  v1 = -23.0;  v2 = 3.255;  v3 = 0.8342;  v4 = -0.2905;  v5 = 0.0;

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