% Script to solve for cloud-topped mixed layer depth, cloudbase and inversion 
% strength given a specified SST, radiative cooling dRN, divergence D, 
% above-inversion theta-v and q profiles.

%  Parameters

  SST = 290;
  D = 5e-6;
  thvp0 = 303;
  dthvpdz = .004;
  qtp0 = .004;
  dqtpdz = 0;
  dRN = 50;
  CTV = 0.01;
  p0 = 102000;
  pref = 97700;
  Tref = SST - 4.5;
  
%  solve for BL depth h using quadratic a2*h^2 + a1*h + a0 = 0.

  q0 = qs(p0, SST);
  thv0 = (100000/p0)^(Rd/Cp) * SST * (1 + delta*q0);
  a2 = D*dthvpdz;
  a1 = D*(thvp0 - thv0);
  a0 = -dRN/(rhoref(SST)*Cp);
  h = (-a1 + sqrt(a1^2 - 4*a0*a2))/(2*a2);
  we = D*h;
  dthvi = a0/we;

%  solve for BL thermo properties

  chi = we/(we+CTV);
  qp = qtp0 + dqtpdz*h;
  qtm = (1 - chi)*q0 + chi*qp;
  thm = thv0/(1 + delta*qtm);
  Tm0 = thm*(p0/1e5)^(Rd/Cp);
  qsm0 = qs(p0, Tm0);
  zb = (qtm - qsm0)/dqsdzu(pref,Tref);


  