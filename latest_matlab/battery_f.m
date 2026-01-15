function xdot = battery_f(x, iT, p)
% x = [soc; vTS; vTL]
soc = x(1);
vTS = x(2);
vTL = x(3);

RS  = p.RS.a  * exp(p.RS.b  * soc) + p.RS.c;
RTS = p.RTS.a * exp(p.RTS.b * soc) + p.RTS.c;
RTL = p.RTL.a * exp(p.RTL.b * soc) + p.RTL.c;

CTS = p.CTS.a * exp(p.CTS.b * soc) + p.CTS.c;
CTL = p.CTL.a * exp(p.CTL.b * soc) + p.CTL.c;

soc_dot = -iT / p.CBat;
vTS_dot = -(vTS / (CTS*RTS)) + iT/CTS;
vTL_dot = -(vTL / (CTL*RTL)) + iT/CTL;

xdot = [soc_dot; vTS_dot; vTL_dot];
end
