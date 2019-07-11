function rho = rhoref(SST)

%  Reference BL density at reference pressure pref and a temp. of SST - 4.5 K

  global pref
  rho =  rhodry(pref, (SST - 4.5));
