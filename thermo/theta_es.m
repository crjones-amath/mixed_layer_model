function thetaes = theta_es(p,Temp)

%  theta_es(p,Temp) is saturation theta-e using thetae and qs formulas from
%   Bolton, 1980, MWR, 108, 1046-1053.

  thetaes = theta_e(p,Temp,qs(p,Temp));