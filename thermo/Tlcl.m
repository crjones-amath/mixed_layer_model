function Tl = Tlcl(p,Temp,qv)

%  Temperature at the LCL
%   From Bolton, 1980, MWR, 108, 1046-1053.

  global Rd Rv
  ev = p.*qv./(Rd/Rv + qv);
  Tl = 2840./(3.5*log(Temp) - log(.01*ev) - 4.805) + 55;
