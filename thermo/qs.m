function qsat = qs(p,Temp)

%  qs(p,Temp) is saturation mixing ratio based on Wexler's formula for es.
%   From Bolton, 1980, MWR, 108, 1046-1053.

%  global Rd Rv
  Rd    = 287.04;
  Rv    = 461.50;

  esat = es(Temp);
  qsat = (Rd/Rv)*esat./(p-esat);
