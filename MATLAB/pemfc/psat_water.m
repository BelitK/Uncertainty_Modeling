function Psat = psat_water(T)
%PSAT_WATER Saturation vapor pressure of water [Pa] at temperature T [K]
% Simple empirical formula valid roughly for 0 to 100 degC

    % Tetens-like formula
    Psat = 610.94 .* exp(17.625 .* (T - 273.15) ./ (T - 30.11));
end
