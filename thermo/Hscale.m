function H = Hscale(Temp)

% Hscale(Temp) = Rd*Temp/g  is pressure scale height

  global Rd g
  H = Rd*Temp/g;