function eta = etag(Temp)

%  etag(Temp) = g/(Cp*Temp)

  global g Cp
  eta = g./(Cp*Temp);
