function qsat = qs_simple(p,Temp)

%  qs_simple(p,Temp) is saturation mixing ratio computed with approximations
%  that L is constant and es << p.

  global L Rv
  qsat = (380./p) .* exp( (L/Rv) .* (1./273. - 1./Temp));
