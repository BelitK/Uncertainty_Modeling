function make_jacobian_abcd(theta, s0, u0, vTS0, vTL0)
% make_jacobian_abcd(theta, s0, u0, vTS0, vTL0)
% Computes Jacobian linearization matrices A,B,C,D and operating output y0
% for x=[soc; vTS; vTL], u=iT, y=vT.


s0 = 0.4;  u0 = 0;  vTS0 = 0;  vTL0 = 0;

% --- unpack theta ---
CBat  = theta(1);

RS_a  = theta(2);  RS_b  = theta(3);  RS_c  = theta(4);

RTS_a = theta(5);  RTS_b = theta(6);  RTS_c = theta(7);
CTS_a = theta(8);  CTS_b = theta(9);  CTS_c = theta(10);

RTL_a = theta(11); RTL_b = theta(12); RTL_c = theta(13);
CTL_a = theta(14); CTL_b = theta(15); CTL_c = theta(16);

v0 = theta(17); v1 = theta(18); v2 = theta(19);
v3 = theta(20); v4 = theta(21); v5 = theta(22);

% --- parameter values at s0 ---
RS  = RS_a  * exp(RS_b  * s0) + RS_c;
RTS = RTS_a * exp(RTS_b * s0) + RTS_c;
RTL = RTL_a * exp(RTL_b * s0) + RTL_c;

CTS = CTS_a * exp(CTS_b * s0) + CTS_c;
CTL = CTL_a * exp(CTL_b * s0) + CTL_c;

% --- derivatives at s0 ---
RS_p  = RS_a  * RS_b  * exp(RS_b  * s0);
RTS_p = RTS_a * RTS_b * exp(RTS_b * s0);
RTL_p = RTL_a * RTL_b * exp(RTL_b * s0);

CTS_p = CTS_a * CTS_b * exp(CTS_b * s0);
CTL_p = CTL_a * CTL_b * exp(CTL_b * s0);

voc   = v0 * exp(v1*s0) + v2 + v3*s0 + v4*s0^2 + v5*s0^3;
voc_p = v0*v1*exp(v1*s0) + v3 + 2*v4*s0 + 3*v5*s0^2;

% --- k terms and derivatives ---
kTS = 1/(CTS*RTS);
kTL = 1/(CTL*RTL);

kTS_p = -(CTS_p*RTS + CTS*RTS_p) / (CTS*RTS)^2;
kTL_p = -(CTL_p*RTL + CTL*RTL_p) / (CTL*RTL)^2;

% --- Build A ---
A = zeros(3,3);

% ds/dt = -u/CBat  -> no x dependence
A(1,:) = [0 0 0];

% dvTS/dt = -vTS*kTS(s) + u*(1/CTS(s))
A(2,1) = -vTS0*kTS_p - u0*(CTS_p/(CTS^2));
A(2,2) = -kTS;
A(2,3) = 0;

% dvTL/dt = -vTL*kTL(s) + u*(1/CTL(s))
A(3,1) = -vTL0*kTL_p - u0*(CTL_p/(CTL^2));
A(3,2) = 0;
A(3,3) = -kTL;

% --- Build B ---
B = [ -1/CBat;
       1/CTS;
       1/CTL ];

% --- Build C, D for y = voc(s) - vTS - vTL - u*RS(s)
C = [ voc_p - u0*RS_p,  -1,  -1 ];
D = -RS;

% --- operating output y0 = h(x0,u0) ---
y0 = voc - vTS0 - vTL0 - u0*RS;

% Export to base workspace for Simulink blocks
assignin('base','A',A);
assignin('base','B',B);
assignin('base','C',C);
assignin('base','D',D);
assignin('base','y0',y0);
assignin('base','x0',[s0; vTS0; vTL0]);
assignin('base','u0',u0);
end
