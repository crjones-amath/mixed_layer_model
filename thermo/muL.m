function m = muL(Temp)

%  muL(Temp) =  (1 - epsi(Temp)*delta)*Lv(Temp)
%   Thermo parameter relating unsaturated sv flux to qt flux.  

  global delta
  m = (1 - epsi(Temp)*delta).*Lv(Temp);
