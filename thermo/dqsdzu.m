function dqsatdzu = dqsdzu(p,Temp)

%  dqsdzu = dqs/dz on a dry adiabat

  global Rd Cp
  dqsatdzu = ((Rd*Temp/Cp) .* dqsdT(p,Temp) + p .* dqsdp(p,Temp))./Hscale(Temp);
