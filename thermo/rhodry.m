function rho = rhodry(p,Temp)

%  rhodry(p,Temp) = p/(Rd*Temp) is density of dry air

  global Rd
  rho = p./(Rd*Temp);
