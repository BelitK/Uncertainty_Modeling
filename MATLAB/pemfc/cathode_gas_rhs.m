function dPc = cathode_gas_rhs(t, Pc, par, in)
%CATHODE_GAS_RHS Reduced cathode gas channel model
%
%   Pc(1) = Pc_O2  [Pa]
%   Pc(2) = Pc_H2O [Pa]
%
%   dPc is a 2x1 vector [dPc_O2; dPc_H2O]

    Pc_O2  = Pc(1);
    Pc_H2O = Pc(2);

    R        = par.R;
    F        = par.F;
    Scell    = par.Scell;
    Ncell    = par.Ncell;
    Vc       = par.Vc;
    Pc_tot   = par.Pc_tot;
    X_O2_air = par.X_O2_air;   % typically 0.21

    Win_c = in.Win_c;   % [mol/s]
    HR_c  = in.HR_c;    % relative humidity
    T     = in.T;       % [K]
    j     = in.j;       % [A/m^2]
    Jw_c  = in.Jw_c;    % [mol/(m^2 s)]

    % Saturation and inlet mole fractions
    Psat      = psat_water(T);            % [Pa]
    Xc_in_H2O = HR_c * Psat / Pc_tot;     % inlet water fraction
    Xc_in_O2  = X_O2_air * (1 - Xc_in_H2O);

    % Outlet molar flow according to reduced model
    % O2 consumption: j*Scell*Ncell / (4F)
    % H2O production: j*Scell*Ncell / (2F)
    % Membrane water flux to cathode: Jw_c*Scell*Ncell
    W_c_out = Win_c ...
        + (Scell * Ncell * j) / (2 * F) ...
        - Scell * Ncell * Jw_c;

    % Oxygen balance
    dPc_O2 = (R * T / Vc) * ( ...
        Win_c * Xc_in_O2 ...
        - W_c_out * (Pc_O2 / Pc_tot) ...
        - (Scell * Ncell * j) / (4 * F) );

    % Water balance
    dPc_H2O = (R * T / Vc) * ( ...
        Win_c * Xc_in_H2O ...
        - W_c_out * (Pc_H2O / Pc_tot) ...
        + (Scell * Ncell * j) / (2 * F) ...
        + Scell * Ncell * Jw_c );

    dPc = [dPc_O2; dPc_H2O];
end
