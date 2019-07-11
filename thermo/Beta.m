function B = Beta(p,Temp)

%  Beta(p,Temp)= (1 + Gamma(p,Temp)*epsi(Temp)*(delta + 1))/(1 + Gamma(p,Temp))
%   Moist thermo parameter relating saturated sv flux to h flux.

  global delta
  B = (1 + Gamma(p,Temp).*epsi(Temp)*(delta + 1))./(1 + Gamma(p,Temp));
