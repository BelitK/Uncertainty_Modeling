function vT = battery_h(x, iT, p)
soc = x(1);
vTS = x(2);
vTL = x(3);

RS = p.RS.a * exp(p.RS.b * soc) + p.RS.c;

voc = p.OCV.v0 * exp(p.OCV.v1 * soc) + ...
      p.OCV.v2 + p.OCV.v3*soc + p.OCV.v4*soc^2 + p.OCV.v5*soc^3;

vT = voc - vTS - vTL - iT * RS;
end
