function dPa = anode_gas_rhs(t, Pa, par, in)
%ANODE_GAS_RHS Reduced anode gas channel model
%
%   Pa(1) = Pa_H2  [Pa]
%   Pa(2) = Pa_H2O [Pa]
%
%   dPa is a 2x1 vector [dPa_H2; dPa_H2O]

    Pa_H2  = Pa(1);
    Pa_H2O = Pa(2);

    R      = par.R;
    F      = par.F;
    Scell  = par.Scell;
    Ncell  = par.Ncell;
    Va     = par.Va;
    Pa_tot = par.Pa_tot;

    Win_a = in.Win_a;  % [mol/s]
    HR_a  = in.HR_a;   % relative humidity
    T     = in.T;      % [K]
    j     = in.j;      % [A/m^2]
    Jw_a  = in.Jw_a;   % [mol/(m^2 s)]

    % Saturation and inlet mole fractions
    Psat      = psat_water(T);            % [Pa]
    Xa_in_H2O = HR_a * Psat / Pa_tot;
    Xa_in_H2  = 1 - Xa_in_H2O;

    % Outlet molar flow according to reduced model
    % H2 consumption: j*Scell*Ncell / (2F)
    % Membrane water flux to anode: Jw_a*Scell*Ncell
    W_a_out = Win_a ...
        - (Scell * Ncell * j) / (2 * F) ...
        - Scell * Ncell * Jw_a;

    % Hydrogen balance
    dPa_H2 = (R * T / Va) * ( ...
        Win_a * Xa_in_H2 ...
        - W_a_out * (Pa_H2 / Pa_tot) ...
        - (Scell * Ncell * j) / (2 * F) );

    % Water balance
    dPa_H2O = (R * T / Va) * ( ...
        Win_a * Xa_in_H2O ...
        - W_a_out * (Pa_H2O / Pa_tot) ...
        + Scell * Ncell * Jw_a );

    dPa = [dPa_H2; dPa_H2O];
end
