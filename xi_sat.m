function [xi, qsat, hf, T] = xi_sat(htrop,hcld,qtrop,qcld,qt,zi,p,varargin)
% calculate mixing fraction (xi) of above-PBL air in an exactly saturated 
% mixture of cloudy + above-PBL air:
global Cp L g

% default options
do_simple = 0; qtol = 1e-7; htol = 1e-4; 

numvarargs = length(varargin);
optargs = {do_simple, qtol, htol};
optargs(1:numvarargs) = varargin;
[do_simple, qtol, htol] = optargs{:};

% first guess for iteration;
qv = qt - qcld; % vapor content at cloud top

xi   = 0.20; 
maxit = 200;

CpTp = (htrop - g.*zi - L.*qtrop);
CpTm = (hcld - g.*zi - L.*qv);

% if htrop < hcld, dxi sign = +1 (changes in same dir as dh)
% else, dxi sign = -1 (changes in opposite dir);
% dxi = 0.25.*sign(hcld - htrop);
dxi = 0.02;

for kk=1:maxit
   CpT  = xi.*(CpTp) + (1-xi).*CpTm - (1-xi).*L.*qcld;
   T    = CpT./Cp;
   if do_simple
      qsat = qs_simple(p,T);
   else
      qsat = qs(p,T);
   end
   hsat = CpT + g.*zi + L.*qsat;
   hf   = xi.*htrop + (1-xi).*hcld; % MSE at saturation
   qxi  = xi.*qtrop + (1-xi).*qt;   % qsat
   dq   = qxi-qsat;
   % use known qsat to find T:
   if (abs(hsat-hf)<htol && abs(qsat-qxi)<qtol)
     % disp(['xi converged after ',num2str(kk),' iterations'])
      break;
   else
      xi = xi+sign(dq).*dxi;
%      dxi = dxi./2;
      dxi = dxi.*(1 - 0.25.*(1+sign(dq))); % halve the step-size when dx increases, not when decreases
   end
   if kk==maxit
      warning('xi_sat:conv','Did not converge')
      disp(['dxi = ', num2str(dxi)])
   end
end