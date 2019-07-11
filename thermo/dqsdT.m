function dqsatdT = dqsdT(p,Temp)

%  dqsdT = d(qs(p, Temp))/d(Temp) = qs(p,Temp)*Lv(Temp)/(Rv*Temp**2)  [K^{-1}]
%   using Wexler and Clausius-Clapeyron formulas

  global Rv
  dqsatdT = qs(p,Temp).*Lv(Temp)./(Rv*Temp.^2);
