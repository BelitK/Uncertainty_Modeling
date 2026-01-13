function [x_next, j] = pemfc_step_paper_clean(x, u)
%#codegen
% Stage 1: solve j from voltage equation only (paper Eq. 1-9 backbone)
% x = [lambda_a; lambda_m; lambda_c]
% u = [V; Win_c; Win_a; Pc; Pa; HRc; HRa; T; Ts]

V  = u(1);
Pc = u(4);
Pa = u(5);
HRc = u(6);
T  = u(8);
Ts = u(9);

% keep x fixed for now
x_next = x;

% constants from paper
R  = 8.314;
F  = 96485.3329;
Pref   = 1.013e5;
Tref   = 298.15;
Eact   = 66000;
Ememb  = 10542;
alpha  = 0.5635;
gamma  = 1.0;
j0ref  = 1.015e-3;
PO2ref = 0.21e5;

dH = -241.83e3;
dS = -163.34;

xm = 25e-6;   % m
xGDL = 15e-6; % m

% O2 diffusion placeholder (must be >0)
eps_gdl = 0.6; tau_gdl = 2.0; DO2 = 2e-5;
DO2eff = (eps_gdl/(tau_gdl*tau_gdl))*DO2;

% inlet partial pressures (starter)
Psat = psat_water_Pa(T);
PcH2O = min(max(HRc*Psat, 1.0), 0.99*Psat);
PcO2  = 0.21*(Pc - PcH2O);     % simple inlet mix, positive if PcH2O reasonable
PcO2  = max(PcO2, 1.0);
PaH2  = 0.9*Pa; PaH2 = max(PaH2, 1.0); % starter

% voltage pieces (paper-like)
Erev = -((dH/(2*F)) - (T*dS/(2*F))) + (R*T/(2*F))*log( ...
       (PaH2/Pref) * (PcO2/Pref)^0.5 * (PcH2O/max(Psat,1.0)) );


disp(['Erev = ', num2str(Erev), ...
          ', PcO2 = ', num2str(PcO2), ...
          ', PcH2O = ', num2str(PcH2O), ...
          ', Psat = ', num2str(Psat)]);

sigma = (0.005139*x(2) - 0.00326) * exp((Ememb/R)*(1/Tref - 1/T));
sigma = max(sigma, 1e-6);
xm_cm = xm*100;

% solve j with a bracket + bisection (monotonic enough)
j_lo = 1e-6;
j_hi = 2.0; % A/cm^2

f_lo = V_of_j(j_lo);
f_hi = V_of_j(j_hi);

% expand bracket if needed
for k=1:10
    if f_lo*f_hi < 0
        break
    end
    j_hi = j_hi*1.5;
    f_hi = V_of_j(j_hi);
end

% if still no bracket, fall back to small j
if f_lo*f_hi > 0
    j = 0.1;
    return
end

% bisection
for it=1:40
    jm = 0.5*(j_lo + j_hi);
    fm = V_of_j(jm);
    if abs(fm) < 1e-6
        j_lo = jm; j_hi = jm;
        break
    end
    if f_lo*fm < 0
        j_hi = jm; f_hi = fm;
    else
        j_lo = jm; f_lo = fm;
    end
end

j = 0.5*(j_lo + j_hi);

    function val = V_of_j(jv)
        % activation
        j0 = j0ref * (PcO2/PO2ref)^gamma * exp((Eact/R)*(1/Tref - 1/T));
        j0 = max(j0, 1e-20);
        eta_act = (R*T/(alpha*F)) * log(max(jv,1e-12)/j0);

        % concentration
        jL_SI = (4*F*DO2eff*PcO2)/(R*T*xGDL); % A/m^2
        jL = max(jL_SI/1e4, 1e-6);           % A/cm^2
        frac = max(1 - jv/jL, 1e-9);
        eta_conc = -(R*T/(alpha*F))*log(frac);

        % ohmic
        eta_ohm = jv * (xm_cm / sigma);

        Vmodel = (Erev - eta_act - eta_conc - eta_ohm);
        val = Vmodel - V; % root when 0
    end

end

function Ps = psat_water_Pa(TK)
TC = TK - 273.15;
TC = min(max(TC, 0), 100);
A = 8.07131; B = 1730.63; C = 233.426;
PmmHg = 10^(A - B/(C + TC));
Ps = max(PmmHg * 133.322, 1.0);
end
