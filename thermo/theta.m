function th = theta(p,Temp)

% Potential temperature Temp*(1e5/p)^(Rd/Cp)

  global Rd Cp
  th = Temp.*(1e5./p).^(Rd/Cp);
