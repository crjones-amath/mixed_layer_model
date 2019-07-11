function [dydt,zb,T0,Tl,zref,qs_surf,hs_surf,dqdzu,dqdzs,rref,ep,mu,B,sv0,svl,...
   ql,sv,dqi,dhi,dsvi,we, W0, Wi, E0, Ei, wstar, BIR, rho,...
   wb_bar, wsv_bar, wh_bar, wq_bar,zgrid,Frad,wsed,driz,dEdz,dWdz,...
   lwp_zb,qv,P, LHF, SHF, p, rhod, ib, iFT, hz, qs0,I0,I1,alpha,lwp,A] = ...
   BW97_MLM_rhs_v2(t,y,QT_TROP,HP,D,SST,SLP,CTh,CTm,zmax,Nd,...
   do_sed,do_full,a1,a2,ased,LES_tuned,Junya,Pref,Tref,do_driz, do_simple, varargin)
%BW97_MLM_rhs    Mixed-layer model "right hand side" function
%
% Mixed Layer Model following
% Bretherton and Wyant, 1997
%
%  Author: Chris Jones
%    Date: 5/28/2013
% Version: 1.1
%-----------------------------------------
% Specify ODEs for
%    qt = ql + qv (total water mixing ratio)
%    h  = Cp T + g z + L qv (moist static energy)
%    zi (inversion height)
% that are suitable to be solved with MATLAB's ODESOLVER functions
%-----------------------------------------
% Input arguments:
%  t : time [s]
%  y : current BL state y = [zi;h;qt].  Units: [m,J/kg,kg/kg]
%  QT_TROP : Free troposphere moisture.  Currently this can either be a
%    constant value or a time-varying function handle. Units: [kg/kg]
%  HP : (function handle) Free troposphere moist static energy profile.
%    May specify HP as a function of z alone, or z and t:  HP(z) or HP(z,t)
%  D : Large scale divergence, ws ~ -D*zi, where ws = subsidence.
%      D is overloaded to allow you to specific a non-linear subsidence
%      function.  If D is a function handle, then ws = D(zi) will be used.
%  SST : Sea surface temperature.  Can be either constant or time-varying
%      function sst = SST(t).  Units: K
%  SLP : Sea-level pressure, assumed to be constant.  Units: Pa
%  CTh,CTm : surface transfer coefficients for h and qt surface fluxes,
%    respectively.  These already include "V" -- i.e., CTh = C*Vsurf, with
%    units [m/s].
%  zmax : Top of grid used for radiation and precipitation.
%  Nd : Droplet concentration [1/m3]
%  do_sed : [1] Include droplet sedimentation in precipitation; 
%           [0] Do not include cloud droplet sedimentation in preciptation
%  do_full: [1] Use full expression for moist therodynamics in BL state
%           [0] Linearized expression for moist thermodynamics in BL state
%  a1,a2,ased:  Parameters in the default TN87-style entrainment closure
%  LES_tuned: [1] Use LES-tuned cloud-base drizzle following Uchida et al
%             [0] Default drizzle, entrainment closure
%  Junya:  [1] Combined with do_full = 0 and do_simple = 1 should reproduce
%              the behavior of the MLM used in Uchida et al (2011).
%          [0] Default behavior.
%  Pref: Reference pressure, suggested to be SLP - 50 hPa.
%  do_driz: [1] Drizzle on. [0] Drizzle off.
%  do_simple: [1] Use qs_simple for saturation; also, use linearized
%    expression for evaporative enhancement in entrainment calculation
%    [0] use qs for saturation; if do_full = 1, iteratively determine
%    evaporative enhancement 'chi^*' value.
%-----------------------------------------
% Optional input arguments (varargin), along with default values:
%
% rad_opts = []; % structure for passing radiation options / other opts
% hadvh = 0; % horizontal advection of h, assumed to be constant w/height
% hadvq = 0; % horizontal advection of qt, assumed to be constant w/height
% hadvzi = 0; % horizontal advection of zi
% do_fixed_LHF = 0; % if do_fixed_LHF = k, then LHF = k W/m2 = constant is
% used.  If do_fixed_LHF = 0, then interactive surface fluxes are used
% do_fixed_SHF = 0; % as in do_fixed_LHF, but for SHF
% debug = 0; % if debug = 1, additional diagnostic output is printed
% we_fun = @calc_we_v2; Function to calculate entrainment. See calc_we_v2
%           for required arguments.
% do_thermo_only = 0; if do_thermo_only = 1, then hold dzi/dt = 0.
%-----------------------------------------
% Output arguments:
%-----------------------------------------
% Calculating the RHS proceeds in the following steps
% 1. Preliminaries: set-up external forcings for given time step
% 2. Calculate atmospheric state, including therodynamic parameters,
%    quanitities, and especially cloud base (zb)
% 3. Call entrainment package to determine we, radiation, and drizzle terms
%    for right hand side source terms.
% 4. Construct RHS.
%-----------------------------------------
global L Cp g gam

%-----------------------------------------
% default options for varargin
%-----------------------------------------
rad_opts = []; % structure for entering radiation options / other opts
hadvh = 0; hadvq = 0; % horizontal advection, assumed for the moment to be constant 
hadvzi = 0;
do_fixed_LHF = 0;
do_fixed_SHF = 0;
debug = 0;
we_fun = @calc_we_v2;
do_thermo_only = 0;

% read in optional arguments
numvarargs = length(varargin);
optargs = {rad_opts,hadvh,hadvq,hadvzi,do_fixed_LHF,...
   do_fixed_SHF,debug,we_fun,do_thermo_only};
optargs(1:numvarargs) = varargin;
[rad_opts,hadvh,hadvq,hadvzi,do_fixed_LHF,do_fixed_SHF,...
   debug,we_fun,do_thermo_only] = optargs{:};

if debug
   disp(['t = ',num2str(t)])
   rad_opts
end

% preliminaries
zi = y(1);
hm = y(2);
qt = y(3);

if isa(SST,'function_handle')
   sst = SST(t);
elseif isa(SST,'numeric')
   sst = SST;
else
   error('BW97_MLM_rhs:SST','SST must be function handle or scalar')
end

if isa(QT_TROP,'function_handle')
   qt_trop = QT_TROP(t);
elseif isa(QT_TROP,'numeric')
   qt_trop = QT_TROP;
else
   error('BW97_MLM_rhs:qt_trop','QT_TROP must be function handle or scalar')
end

if nargin(HP)==2
   hp = @(z) HP(z,t);
elseif nargin(HP)==1
   hp = @(z) HP(z);
else
   warning('BW97_MLM_rhs2:hp','argument problem with hp')
end

% overload D, such that if D is a function it refers to ws(zi).
% if D = constant, then ws(zi) = D*zi
% CONVECTION: D = POSITIVE or ws(zi) = NEGATIVE
if isa(D,'function_handle')
    ws = D(zi); %D(zi) should be a negative number
else
    ws = -D.*zi; % D should be a positive number
end

Tref = sst - 5;

% calculate boundary layer state
%  - calculate zi
%  - construction zgrid.
%  - Indices: zgrid(ib) = zb; zgrid(iFT-1) = zi;
if do_full
   do_double_zi = 0; % 0 => use an offset between zi- and zi+
   [p, T, qv, qsat, ql, Tv, rho, zb, zgrid, qs_surf, hs_surf,dqdzu,...
   dqdzs,T0, Tl, qs0, Theta, Thetav, s, sv, sl,rref,ri,ib,iFT,...
   mu, B, ep, sig, hz, rhod, sv0, qtz, zref,svl] = ...
   MLM_BL_state_full(hm,qt,zi,hp,qt_trop,sst,SLP,zmax,Tref,Pref,Junya,...
      do_simple, do_double_zi);

   rho_e = rhod; % use rho_dry in we calc, which plays into lwp and wsed
   qs0   = qsat(1);
else
   do_double_zi = 0;
   [T0, qs0, qs_surf, hs_surf, rref, dqdzu, dqdzs,...
   zb, ib, iFT, qtz, ql, qv, hz, T, Tv, ri, p, qsat, Theta, Thetav,...
   s, sv, sl, svl, mu, B, ep, sig, sv0, zgrid, zref, Tl] = ...
   MLM_BL_state_linearized(zi,qt,hm,qt_trop,hp,SLP,sst,zmax,Tref,Pref,Junya,do_simple,do_double_zi);

   rho   = rref; rhod = rref;
   rho_e = rref; % use fixed rref in we calc
end

%---------------------
% inversion jumps
%---------------------
% dqi  = qt_trop - qt;
% dhi  = hp(zi)  - hm;
% dsvi = svp(zi) - sv(zi);

dqi  = qtz(iFT) - qtz(iFT-1);
dhi  = hz(iFT) - hz(iFT-1);
dsvi = sv(iFT) - sv(iFT-1);

%---------------------------
% radiation options
% rad_opts.t0 specifies start time to apply coefficients in rad_opts
%---------------------------
rad_we = rad_opts;
if ~isempty(rad_opts)
   if t<=rad_opts.t0
      rad_we.F0 = [];
      rad_we.F1 = [];
      rad_we.k  = [];
      rad_we.coef3 = [];
      if isfield(rad_we,'wefixed')
         rad_we.wefixed = rad_opts.wefixed(1);
      end
      if isfield(rad_we,'we_propto_zc')
         rad_we.we_propto_zc = rad_opts.we_propto_zc(1);
      end
   else
      if isfield(rad_we,'wefixed')
         rad_we.wefixed = rad_opts.wefixed(2);
      end
      if isfield(rad_we,'we_propto_zc')
         rad_we.we_propto_zc = rad_opts.we_propto_zc(2);
      end      
   end
   % specified bl cool
   if isfield(rad_we,'hadvs_Wm2')
       hadvq = 0;
       hadvh = rad_we.hadvs_Wm2/(zi*rref);
   end
end

% entrainment, radiation, and precipitation calculations:
[we, W0, Wi, E0, Ei, wstar, BIR,...
   wb_bar, wsv_bar, wh_bar, wq_bar,Frad,wsed,driz,lwp_zb,P,LHF,SHF,I0,I1,alpha,lwp,A] = ...
   we_fun(zb,zi,zgrid, ql, rho_e, rref, Nd, sv0, dsvi, hm, qt, dhi, dqi,...
   hs_surf, qs_surf, mu, B, ep, CTh, CTm, D, do_sed, ib, iFT, ...
   a1, a2, ased,p,LES_tuned,do_full,Junya,do_driz,do_simple, rad_we,...
   do_fixed_LHF,do_fixed_SHF,debug);

dEdz = (Ei - E0)./zi; % eq. (4)
dWdz = (Wi-W0)./zi; % eq. (5)

% finally we have the RHS:
dzi = we + ws   - hadvzi; % eq (1)
dhm =    - dEdz - hadvh; % eq (2)
dqt =    - dWdz - hadvq; % eq (3)

dydt = [dzi; dhm; dqt];

if do_thermo_only
    dydt = [0; dhm; dqt];
end

% debug purposes:
if debug
    disp(['t = ',num2str(t/3600),' h'])
    disp(['zi = ',num2str(zi)])
    disp(['zb = ',num2str(zb)])
    disp(['BIR = ',num2str(BIR)])
    disp(['we = ',num2str(we)])
end

end