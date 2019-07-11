function thetae = theta_e(p,Temp,qv)

%  theta_e(p,Temp,qv) from eqn. 43 of Bolton, 1980, MWR, 108, 1046-1053.

   global Rd Cp 
   Tl = Tlcl(p,Temp,qv);
   Pi = (p/100000).^((Rd/Cp).*(1-0.28*qv));
   thetae = Temp./Pi .* exp((3376./Tl - 2.54).*qv.*(1 + 0.81*qv));
