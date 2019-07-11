function dqsatdp = dqsdp(p,Temp)

%  dqsdp(p,Temp) = d(qs)/d(Temp) = -qs(p,Temp)/(p - es(Temp))   [Pa^{-1}]

  dqsatdp = -qs(p,Temp)./(p - es(Temp));
