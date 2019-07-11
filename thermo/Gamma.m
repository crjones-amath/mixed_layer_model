function gam = Gamma(p,Temp)

%  Gamma(p,Temp) = (Lv(Temp)/Cp)*dqsdT(p,Temp)

  global Cp
  gam =  (Lv(Temp)/Cp).*dqsdT(p,Temp);
