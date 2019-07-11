function epsilon = epsi(Temp)

%  epsi(Temp) = Cp*Temp/Lv(Temp)

  global Cp
  epsilon = Cp*Temp./Lv(Temp);
