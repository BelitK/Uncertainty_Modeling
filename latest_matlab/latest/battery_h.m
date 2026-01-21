function vT = battery_h(x, i, theta)
sigma = x(1);
vTS   = x(2);
vTL   = x(3);

sigma = min(max(sigma,0),1);

v0 = theta(2); v1 = theta(3); v2 = theta(4);
v3 = theta(5); v4 = theta(6); v5 = theta(7);

RS_a = theta(8); RS_b = theta(9); RS_c = theta(10);

vOC = v0*exp(v1*sigma) + v2 + v3*sigma + v4*sigma^2 + v5*sigma^3;
RS  = RS_a*exp(RS_b*sigma) + RS_c;
RS  = max(RS, 1e-12);

vT = vOC - vTS - vTL - i*RS;
end
