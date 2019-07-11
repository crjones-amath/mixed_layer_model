function [p, T, qv, qsat, ql, Tv, rho, zb, z, qs_surf, hs_surf,dqdzu,...
      dqdzs,T0, Tl, qs0, Theta, Thetav, s, sv, sl,rref,ri,ib,iFT,...
      mu, B, eps, sig, hz, rhod, sv0, qtz, zref, svl]...
   = MLM_BL_state_full(h,qt,zi,hp,qtp,TS,PS,zmax, Tref, Pref,Junya,...
      do_simple, do_double_zi)
% Calculates full diagnostic state of BL for mixed layer
% using either linearized or fully nonlinear moist thermodynamics
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
global g Rd Rv Cp L delta

if do_simple
   qsf = @qs_simple; % use simple q*
else
   qsf = @qs; % us full q*
end

% Surface properties:
T0   = (h - L.*qt)./Cp;   % air temp at surface, assumes qt < qs at surface
qs0  = qsf(PS,T0);         % saturation q at surface (kg/kg)
qs_surf = qsf(PS,TS);        % saturation of q at sea surface (kg/kg)

hs_surf = Cp.*TS+L.*qs_surf;   % MSE at SLP, SST, and qsat

% Tref = TS-4.5;
% Pref = PS-45.*100; % ref pressure
rref = rhodry(Pref,Tref);

zref = (h - Cp.*Tref -L.*qt) ./ g; % needed to see if zref > zb or not

[p, T, qv, qsat, ql, Tv, rho, zb, z, ib, iFT, qtz, rhod,hz] = ...
   BL_profile_full(PS,qt,h,zi,zmax,qtp,hp,qsf, do_double_zi,rref);
% [p, T, qv, qsat, ql, Tv, rho, zb, z,ib,iFT] = BL_profile(PS,qt,h,zi,zmax,qtp);
ri  = rho(iFT-1); % reference density just below inversion

dqdzu = dqsdzu(p,T); % qsat lapse rate along moist adiabat
dqdzs = dqsdzs(p,T); % qsat lapse rate along moist adiabat
Tl = Tlcl(PS,T0,qt); % LCL temperature

% other thermodynamic variables that may be of interest
Theta  = theta(p,T);
Thetav = theta(p,Tv);
s      = Cp.*T  + g.*z;
sv     = Cp.*Tv + g.*z;
sl     = Cp.*T  + g.*z - L.*ql;
svl    = hz - muL(T).*qtz;

% at reference p,T
% Reference level state variables
mu = muL(Tref)./L;
B  = Beta(Pref,Tref);
eps = epsi(Tref);           % epsilon at reference Temp/pres
sig = B.*mu-eps;

% sv0  = interp1(z(1:iFT-1),sv(1:iFT-1),zref); % gets around double-lev at iFT
sv0  = (h-L.*qt).*(1+0.608.*qt); % This is right?

if Junya
  sv0  = (h-Cp.*qt).*(1+0.608.*qt); % is this right?
end

end

function [p, T, qv, qsat, ql, Tv, rho, zb, z, ib, iFT, qtz, rhod, hz] = ...
   BL_profile_full(p0,qt,h,zi,zmax,qtp,hp,qsf, do_double_zi,rref)
global Cp Rd Rv g L delta

ep    = Rd./Rv;
alpha = Cp./(Rd.*(1+delta.*qt));

%----------------------------------------------------------------
% +++cjones
% edit:  3/6/2012:  Improved -- zb can be determined exactly
%     assuming h,qv=qt constant with height up to zb
%----------------------------------------------------------------
% step 1:  determine zb
sl     = h - L.*qt;
zb_fun = @(z) qsf(p0.*((sl-g.*z)./sl).^alpha,...
   1./Cp.*(sl-g.*z)) - qt;
zb     = max(fzero(zb_fun,zi-150), 0); % zb >= 0.

% disp(['debug:  zi = ',num2str(zi)])
% disp(['debug:  zb = ',num2str(zb)])

% step 2:  setup grid
[zgrid, dum1,dum2,dum3, z] = MLM_grid(zb,zi,zmax,100);
ix = find(z==zi);  % should have length 2
ib  = find(z==zb); % should be single-valued
iFT = ix(end);     % index where FT starts

if ~do_double_zi % shift zi(+) up a bit
   offset = 0.1; % fraction of grid spacing to shift it up
   z(iFT) = z(iFT) + offset.*(z(iFT+1)-z(iFT));
end
zSC = z<=zb;
zC  = z>zb & z<=zi;
zFT = z>zi;

% step 3:  fill in qv, ql, hz, qtz, p, T, and rho:
% Subcloud
qv(zSC) = qt;
ql(zSC) = 0;
T(zSC)  = (h - L.*qv(zSC) - g.*z(zSC))./Cp;
T0 = T(1);
% Tv(zSC) = T(zSC).*(qt+ep)./(ep.*(1+qt));
Tv(zSC) = T(zSC).*(1+delta.*qt);   % this approximation is fine
p(zSC)  = p0.*(T(zSC)./T0).^alpha;

% Free trop
if isa(qtp,'function_handle')
   qv(zFT) = qtp(z(zFT));
else
   qv(zFT) = qtp;
end
ql(zFT) = 0;
T(zFT)  = (hp(z(zFT))-L.*qv(zFT) - g.*z(zFT))./Cp;
Tv(zFT) = T(zFT).*(qv(zFT)+ep)./(ep.*(1+qv(zFT)));
dz = diff(z);

for kk=2:length(z)
   if zC(kk)
      qv(kk) = qv(kk-1)-dqsdzs(p(kk-1),T(kk-1)).*dz(kk-1); % will verify in a sec
      ql(kk) = qt - qv(kk);
      T(kk)  = (h - L.*qv(kk) - g.*z(kk))./Cp;
      Tv(kk) = T(kk).*(qv(kk)+ep)./(ep.*(1+ql(kk)+qv(kk)));
   end
   if (z(kk)==zi && z(kk-1)==zi) %zi, but w/free trop values
      if isa(qtp,'function_handle')
         qv(kk) = qtp(zi);
      else
         qv(kk) = qtp;
      end
      ql(kk) = 0;
      T(kk) = (hp(zi)-L.*qv(kk) - g.*zi)./Cp;
      Tv(kk) = T(kk).*(qv(kk)+ep)./(ep.*(1+ql(kk)+qv(kk)));
   end
   p(kk) = p(kk-1).*exp(-g.*dz(kk-1)./(Rd.*Tv(kk-1))); % slightly more accurate, but shouldn't make a diff
end
   rho = p./(Rd.*Tv);
   qsat = qsf(p,T);
   qtz = qv + ql; % qt on grid
   rhod = rho./(1+qtz); % density of dry air
   
   hz = hp(z);
   hz(1:iFT-1) = h;
end