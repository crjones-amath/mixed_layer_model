function F = rad_simple(z,qc,rho,zi,varargin)
global Cp
%-----------------------------------------
% idealized radiation scheme
% use same as dycoms II RF01 idealized case
% (dF/dz)_rad = k*rho*q_c*(F0 exp(-k*lwp(z)) - F1 exp(-k [LWP_b - lwp(z)]))
% Parameters from Stevens et al (2005)
%-----------------------------------------
% z   = vertical grid
% qc  = cloud liquid water at each z grid point (size = size(z))
% rho = density of dry air
% optional inputs:
%    varargin(1) = D  ( 0 1/s   default)
%    varargin(2) = F0 (70 W/m2  default)
%    varargin(3) = F1 (22 W/m2  default)
%    varargin(4) = k  (85 m2/kg default)
%-----------------------------------------
D  =  0; % By default, last term isn't included
F0 = 70; % W/m2, downwelling radiation
F1 = 22; % W/m2, cloudbase upwelling radiation
k  = 85; % m2/kg
numvarargs = length(varargin);
if numvarargs>4
   error('rad_simple:TooManyArgs','Expect at most 4 optional arguments')
end

% set defaults for optional inputs
default_opts = {D F0 F1 k};
optargs = default_opts;
optargs(1:numvarargs) = varargin;
empty_inputs = cellfun(@isempty,optargs); % replace any empty inputs with defaults
optargs(empty_inputs) = default_opts(empty_inputs); 
[D, F0, F1, k] = optargs{:};

Q  = cumtrapz(z,rho.*qc);
Qc = LWP(z,qc,rho);

F    = F0.*exp(-k.*Qc) + F1.*exp(-k.*Q);
% dFdz = k.*rho.*qc.*(F0.*exp(-k.*lwp) - F1.*exp(-k.*(lwp_zb - lwp)));

%--------------------------------------
% 3rd term in Stevens etal (2005)
%--------------------------------------
if length(rho)==length(z)
   rho_i = interp1(z,rho,zi-1.0); % rho_i = rho just below inversion
else
   rho_i = rho; % use rref until I figure out how to deal with this
end
F = F + (rho_i.*Cp.*D.*( 0.25*(z-zi).^(4/3) + zi.*(z-zi).^(1/3) )).*(z>zi);
end