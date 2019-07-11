function Lvap = Lv(Temp)

%  Lv(Temp) = 2.501e6 + (Cpv-Cw)*(Temp-273)  [J/kg]
%   Latent heat of vaporization with temperature correction 
%   From Bolton, 1980, MWR, 108, 1046-1053.

  global Cpv Cw
  Lvap = 2.501e6 + (Cpv-Cw)*(Temp-273);