function [wsed, P] = sedi_flux(ql,Nd,varargin)
% sedimentation flux for MLM following Uchida, Bretherton, and Blossey (2010)
% wsed  = c (3 / [4 pi rho_ell Nd])^(2/3) * (rho qc)^(5/3) exp(5 ln^2 sig_g)
% ql    = cloud top LWC (kg/kg)
% Nd    = droplet concentration (1/m^3);
%--------------------------------------
% varargin order = rho_a, rho_w, debug
% rho_a = air density (kg/m3)
% rho_w = liquid water density (kg/m3)
% debug = 1 (true), 0 (false)
numvarargs = length(varargin);
   rho_a = 1.2;
   rho_w = 1000;
   debug = 0;
if numvarargs>0 % use full formula
   optargs = {rho_a, rho_w, debug};
   optargs(1:numvarargs) = varargin;
   [rho_a, rho_w, debug] = optargs{:};

   c   = 1.19e8; %(ms)^(-1)
   sig = 1.2; 
   r   = ( (3.*rho_a.*ql)./(4.*pi.*rho_w.*Nd) ).^(1/3);
   wsed = c.*r.^2.*exp(5.*( log(sig).^2)  );
   P    = wsed.*(rho_a.*ql);  % sedimentation flux in kg / m2 s
   
   if debug
   subplot(2,1,1)
   plot(ql(ql>0),r(ql>0))
   ylabel('r (m)')

   subplot(2,1,2)
   plot(ql(ql>0),wsed(ql>0))
   xlabel('q_l')
   ylabel('w_{sed} (m/s)')
   end
else
   wsed = 6e5.*(ql./Nd).^(2/3); 
   P    = wsed.*(rho_a).*ql;
end
end