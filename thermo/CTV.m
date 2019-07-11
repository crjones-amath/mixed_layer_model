function C = CTV(u10,v10)

%  CTV(u10,v10) = (1.e-3)*(1 + 0.07*vtot)* vtot,  
%   vtot = sqrt(u10^2 + v10^2)  [m/s]
%   Ocean surface bulk transfer coefficient from Schubert et al. 1979.

  vtot = sqrt(u10^2 + v10^2);
  C = (1.e-3)*(1 + 0.07*vtot)* vtot;

