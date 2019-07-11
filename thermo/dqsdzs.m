function dqsatdzs = dqsdzs(p,Temp)

%  dqsdzs(p,Temp) = dqs/dz on a moist adiabat.

  dqsatdzs = dqsdzu(p,Temp)./(1 + Gamma(p,Temp));
