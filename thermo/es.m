function esat = es(Temp,flag)

%  es(Temp) = 611.2*exp(17.67*(Temp-273)/(Temp-29.5))  [Pa]
%   Saturation vapor pressure (Wexler's formula)
%   From Bolton, 1980, MWR, 108, 1046-1053.

if ~exist('flag','var')
    flag = 'liq';
end

switch flag
    case 'liq'
        esat = 611.2*exp(17.67*(Temp-273)./(Temp-29.5));
    case 'ice'
        esat = NaN; % fill in soon
    otherwise
        error('es:iceLiqFlag','Specify either ice or liq')
end
