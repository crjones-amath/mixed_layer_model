function Fp = drizzle_flux(z,zb,zi,lwp,Nd,varargin)
% precip flux for MLM following Caldwell/Bretherton 2009
% expected units:
% Fp:  mm/s = kg /(m2 s); rbar: micrometers; Nd:  (1/m^3); LWP:  kg/m^2;
% input:
% z   = z values for Fp (size(Fp) = size(z));
% zb  = cloud base
% zi  = inversion
% lwp = liq water path at cloud base (kg/m^2);
% Nd  = droplet concentration (1./cm^3);
% varargin = LES_tuned, k
%-----------------------------------------------
% unit conversions
mm_per_sec  = 1e3;              % (m/s)   -> (mm/s)
g_per_m2    = 1e3;              % (kg/m2) -> (g/m2)
cm3_per_m3  = 1e6;              % (cm^3)  -> (m^3);

numvarargs = length(varargin);
% default values
LES_tuned = 0;
k         = 320; % Junya's version uses k = 32.0
Junya     = 0;
rbar      = 60;
optargs = {LES_tuned, k, Junya,rbar};
optargs(1:numvarargs) = varargin;
[LES_tuned, k, Junya, rbar] = optargs{:};

lwp = lwp.*g_per_m2;        % convert kg/m2 to g/m2
Nd  = Nd./cm3_per_m3;       % convert 1/m3 to 1/cm3

Fb = -24.*0.0156./86400.*(lwp./Nd).^(1.75);

if LES_tuned
   Fb =-0.023.*((lwp./Nd)^3.25)./86400; %new Chris' param
end

Fp = Fb.*( exp(-k.*( (zb-z)./(rbar.^(2.5)) ).^(1.5)).*(z<zb) ...
   + (1 - ( (z-zb)./(zi-zb) ).^3).*(z>=zb & z <zi) );

if Junya
% Junya's version?
Fp = Fb.*( exp(-k.*( (zb-z)./(rbar.^2) ).^(1.5)).*(z<zb) ...
   + (1 - ( (z-zb)./(zi-zb) ).^3).*(z>=zb & z <zi) );
end