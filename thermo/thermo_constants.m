%  Script defining important global constants (Emanuel, Appendix 2)

  global g stefan r_earth Omega Rd Rv Cp Cpv Cw L delta psurf pref

  g     = 9.8;
  stefan= 5.67e-8;
  r_earth = 6370*1000;
  Omega   = 2*pi/86400;

  Rd    = 287.04;
  Rv    = 461.50;
  Cp    = 1005.7;
  Cpv   = 1870;
  Cw    = 4190;
  L     = 2.5e6;
  delta = 0.608;

%  Reference values and definitions for subtropical CBL's

  psurf = 1022*100;
  pref  = psurf - 45*100;
