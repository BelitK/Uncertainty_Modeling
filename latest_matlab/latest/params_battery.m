% Parameters from "Battery State Observation and Condition Monitoring..."
% Table I: Identified parameters of a new NCR18650A battery
Ts = 0.1;
p = struct();

% Rated capacity: 3100 mAh = 3.1 Ah = 3.1*3600 Coulombs
p.CBat = 3.1 * 3600;

% OCV: voc = v0*exp(v1*sigma) + v2 + v3*sigma + v4*sigma^2 + v5*sigma^3
p.v0 = -1;
p.v1 = -23;
p.v2 = 3.255;
p.v3 = 0.8342;
p.v4 = -0.2905;   % (Table formatting in PDF is a bit messy; this is the listed coefficient)
p.v5 = 0;         % Table I does not show an extra value for v5 in the extracted text

% Series resistor: RS = RSa*exp(RSb*sigma) + RSc
p.RS_a = 0.25;
p.RS_b = -20;
p.RS_c = 0.07;

% Short time constant branch
% RTS = RTSa*exp(RTSb*sigma) + RTSc
p.RTS_a = 1;
p.RTS_b = -30;
p.RTS_c = 0.015;

% CTS = CTSa*exp(CTSb*sigma) + CTSc
p.CTS_a = -900;
p.CTS_b = -2;
p.CTS_c = 1000;

% Long time constant branch
% RTL = RTLa*exp(RTLb*sigma) + RTLc
p.RTL_a = 0.01;
p.RTL_b = -4;
p.RTL_c = 0.05;

% CTL = CTLa*exp(CTLb*sigma) + CTLc
p.CTL_a = 25000;
p.CTL_b = -2;
p.CTL_c = 2000;



theta = zeros(25,1);

theta(1)  = p.CBat;

% OCV (6)
theta(2)  = p.v0;
theta(3)  = p.v1;
theta(4)  = p.v2;
theta(5)  = p.v3;
theta(6)  = p.v4;
theta(7)  = p.v5;

% RS (3)
theta(8)  = p.RS_a;
theta(9)  = p.RS_b;
theta(10) = p.RS_c;

% RTS (3)
theta(11) = p.RTS_a;
theta(12) = p.RTS_b;
theta(13) = p.RTS_c;

% CTS (3)
theta(14) = p.CTS_a;
theta(15) = p.CTS_b;
theta(16) = p.CTS_c;

% RTL (3)
theta(17) = p.RTL_a;
theta(18) = p.RTL_b;
theta(19) = p.RTL_c;

% CTL (3)
theta(20) = p.CTL_a;
theta(21) = p.CTL_b;
theta(22) = p.CTL_c;

% (Optional padding for future if you need it)
% keep remaining as zeros unless you need more params
theta(23) = 0;
theta(24) = 0;
theta(25) = 0;