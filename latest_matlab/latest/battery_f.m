function dx = battery_f(x, i, theta)
sigma = x(1);
vTS   = x(2);
vTL   = x(3);

sigma = min(max(sigma,0),1);

CBat = theta(1);

RTS_a = theta(11); RTS_b = theta(12); RTS_c = theta(13);
CTS_a = theta(14); CTS_b = theta(15); CTS_c = theta(16);
RTL_a = theta(17); RTL_b = theta(18); RTL_c = theta(19);
CTL_a = theta(20); CTL_b = theta(21); CTL_c = theta(22);

RTS = RTS_a*exp(RTS_b*sigma) + RTS_c;
CTS = CTS_a*exp(CTS_b*sigma) + CTS_c;
RTL = RTL_a*exp(RTL_b*sigma) + RTL_c;
CTL = CTL_a*exp(CTL_b*sigma) + CTL_c;

% safety floors
tiny = 1e-12;
RTS = max(RTS, tiny); CTS = max(CTS, tiny);
RTL = max(RTL, tiny); CTL = max(CTL, tiny);

dsigma = -i / CBat;
dvTS   = -(1/(RTS*CTS))*vTS + (1/CTS)*i;
dvTL   = -(1/(RTL*CTL))*vTL + (1/CTL)*i;

dx = [dsigma; dvTS; dvTL];
end
