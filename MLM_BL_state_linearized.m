function [T0, qs0, qs_surf, hs_surf, rref, dqdzu, dqdzs,...
   zb, ib, iFT, qtz, ql, qv, hz, T, Tv, ri, p, qsat, Theta, Thetav,...
   s, sv, sl, svl, mu, B, eps, sig, sv0, z, zref, Tl] = ...
   MLM_BL_state_linearized(zi,qt,h,qtp,hp,PS,TS,zmax, Tref, Pref, ...
   Junya, do_simple, do_double_zi)
% Calculates diagnostic state of BL for mixed layer
% using linearized moist thermodynamics
%------------------------------------------------------------
% Inputs:
% h  = CpT + gz + Lqv (Moist static energy) in mixed layer
% qt = qv + ql        (total water mixing ratio) in mixed layer
% zi = inversion height
% hp = free troposphere MSE profile (function handle)
% qtp = free troposhere qt profile
% SST = sea surface temp (in K)
% SLP = sea level pressure
%------------------------------------------------------------
% Outputs:
%------------------------------------------------------------
% Comments:
% 1. zgrid is same as used for radiation and LWP in MLM
% 2. currently using Cp = const, but could (should?) use
%    Cp = Cp.*(1 + r), where r = qv/(1+qv) = spec. humidity
%------------------------------------------------------------
global g Rd Rv Cp L delta gam
ep = Rd./Rv;

% Surface properties:
T0   = (h - L.*qt)./Cp;   % air temp at surface, assumes qt < qs at surface
qs0  = qs(PS,T0);         % saturation q at surface (kg/kg)
qs_surf = qs(PS,TS);        % saturation of q at sea surface (kg/kg)

gam = Gamma(Pref,Tref);

if do_simple
   qs0     = qs_simple(PS,T0);
   qs_surf = qs_simple(PS,TS);
   qs_ref  = qs_simple(Pref,Tref);
end
hs_surf = Cp.*TS+L.*qs_surf;   % MSE at SLP, SST, and qsat

%-----------------------
% test -- temporary!
% qs0 = qs_surf;

% Reference pressure, temp, and density
% rfact = 5; %4.5 by default, 5 in Junya's code
% Tref = TS-rfact;
% Pref = PS-10.*rfact.*100; % ref pressure
rref = rhodry(Pref,Tref);

% lapse rates
if ~do_simple
   dqdzu = dqsdzu(Pref,Tref); % reference dqs/dz along dry adiabat
   dqdzs = dqsdzs(Pref,Tref); % reference dqs/dz along moist adiabat
else
   dqsdT = qs_ref.*L./(Rv.*Tref.^2);
   dqsdp = -qs_ref./Pref;
   gam   = (L./Cp).*dqsdT;
   H_scale = Rd.*Tref./g;
%   ep2   = Cp.*Tref./L;
%   B     = (1+gam.*ep2.*(delta+1))./(1+gam);
   dqdzu = (  (Rd.*Tref./Cp) .* dqsdT + Pref.*dqsdp)./H_scale;
   dqdzs = dqdzu./(1+gam);
end

dqdz2 = 0.000004;
   
% determine zb
zb    = max(0, (qs0-qt)./dqdzu); % zb must be nonnegative
% zb    = max(0, (qs0-qt)./dqdz2); % zb must be nonnegative

if do_simple
   %calculate zb as Junya does
   zb = ((1+gam).*(qs_surf-qt)-gam.*(hs_surf-h)./L)/(dqdzu);
end

Tl     = (h - g.*zb - L.*qt)./Cp;

% fill in grid
[zgrid, dum1,dum2,dum3, z] = MLM_grid(zb,zi,zmax,100);
ib  = find(z==zb); % zb is a grid point, so this should exist
iFT = find(z==zi,1,'last'); % zi+, where free trop begins

if ~do_double_zi % shift zi(+) up a bit
   offset = 0.1; % fraction of grid spacing to shift it up
   z(iFT) = z(iFT) + offset.*(z(iFT+1)-z(iFT));
end

% fill values of qt, ql, qv, and h on grid, 
% also p, T, Tv, and qsat
   qtz = zeros(size(z));
   qtz(1:iFT-1) = qt; qtz(iFT:end) = qtp;
   ql = 0.*qtz;
   ql(ib:iFT-1) = dqdzs.*(z(ib:iFT-1)-zb); % increases linearly in cloud
   qv = qtz-ql;
   hz = hp(z);
   hz(1:iFT-1) = h;
   
   T    = (hz - g.*z - L.*qv)./Cp;
   p    = PS - rref.*g.*z;
   ri   = rhodry(p(iFT-1),T(iFT-1)); % density near inversion
   p(iFT:end) = p(iFT) - ri.*g.*(z(iFT:end)-zi); % better estimate for free trop?  still never really used
   qsat = 0.*hz; 
   qsat(1:ib)     = qs0-dqdzu.*z(1:ib);
   qsat(ib:iFT-1) = qv(ib:iFT-1);
   qsat(iFT:end)  = qs(p(iFT:end),T(iFT:end)); % never really used
   if do_simple
      qsat(iFT:end) = qs_simple(p(iFT:end),T(iFT:end));
   end
   
% Reference level state variables
mu = muL(Tref)./L;
B  = Beta(Pref,Tref);
eps = epsi(Tref);           % epsilon at reference Temp/pres
sig = B.*mu-eps;

Tv   = T + Tref.*(delta.*qv-ql); % linearized virtual temp
sv   = Cp.*Tv + g.*z;
svl  = hz - mu.*L.*qt;

zref = (h - Cp.*Tref -L.*qt) ./ g; % needed to see if zref > zb or not

%sv0  = interp1(z(1:iFT-1),sv(1:iFT-1),zref);
sv0  = (h-L.*qt).*(1+0.608.*qt); % is this right?!

if Junya
%  sv0  = (h-Cp.*qt).*(1+0.608.*qt); % is this right?!
   sv0  = (h-L.*qt).*(1+0.608.*qt); % is this right?!
end
   
% finally, moist-adiabatically conserved quantities:
% other thermodynamic variables that may be of interest
Theta  = theta(p,T);
Thetav = theta(p,Tv);
s      = Cp.*T  + g.*z;
% sv     = Cp.*Tv + g.*z;
sl     = Cp.*T  + g.*z - L.*ql;

end