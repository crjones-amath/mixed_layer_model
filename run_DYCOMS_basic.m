%% DYCOMSII-RF01 MLM
% MLM follows Bretherton and Wyant 1997
% Simple radiation (Stevens etal 2005)
% Sedimentation follows Uchida et al (2010)
% Drizzle follows Caldwell and Bretherton (2009)

addpath('./thermo/');
load_thermo_and_units; % script to load thermo constants and unit conversions

%% Basics
% See BW97_MLM_rhs_v2 for further description of these options

tf = 5.*days_to_sec; % end time in seconds
do_sed = 1; % do sedimentation
do_driz = 1; % do drizzle
do_full = 0; % (0) linearized thermo, (1) nonlinearized thermodynamics
do_simple = 0; % (1) qs = qs_simple; Schubert-style cloud base calculation
LES_tuned = 0; % (1) Use LES_tuned drizzle, entraiment coefs.
Junya = 0; % set Junya = 1 to be consistent with Uchida et al (2010) MLM
do_hadv = 0; % no horizontal advection
hadvT = 0; % only relevant if do_hadv = 1;

% entrainment closure parameters
a1 = 0.2;
a2 = 25;
ased = 9; 
   
%% Additional options that can be changed
col = 'b';   % color for making plots
alt_y0 = []; % alternate initial conditions
do_fixed_LHF = 0;
do_fixed_SHF = 0;
Ndin = 150; % in 1/cm3
we_funhandle = @calc_we_v2;
rad_opts = [];

%% Main code block
do_lin = ~do_full;  % linear thermodynamics (true) or nonlinear (false)

% Time domain
% overload tf to allow specification of tspan
if length(tf)>1
   tspan = tf;
else
   tspan = 0:3600:tf;      % hourly output
end

zmax  = 3000;           % domain height (m)

% Initial Conditions
qtm     = 9.*gkg;       % mixed layer mixing ratio
zi      = 840;          % inversion height (m)
hp      = @(z) 303.92.*kJkg + 6.*z; % hp(0) + (6 kJ/km) z;

qt_trop = 1.5.*gkg;       % free tropospheric qt (g/kg)
D   = 3.75e-6;            % divergence (1/s)
V   = sqrt(6^2+(4.25)^2); % 10-m wind (m/s)
CTh = 0.001.*V;           % Transfer coefficient for bulk aero. SHF
CTm = 0.001.*V;           % Transfer coeff. for bulk aero. LHF

if ~isempty(rad_opts)
   if isfield(rad_opts,'div')
      D = rad_opts.div;
   end
   if isfield(rad_opts,'hp0')
      hp = @(z) rad_opts.hp0 + 6.*z;
   end
end

% CTm, CTh overloaded to be LHF, SHF when fixed fluxes are used
if do_fixed_LHF
   CTm = 115; % LHF, in W/m2
end
if do_fixed_SHF
   CTh = 15; % SHF, in W/m2
end

SST0 = 292.5;             % Initial SST (K)
SST  = SST0;              % could be a function if SST varies with time
SLP = 1017.8.*hPa;        % Sea level pressure
Nd  = Ndin./cm3;          % droplet concentration (1/m^3);

Pref = SLP-50.*hPa; % reference pressure (Pa)
Tref = SST0 - 5;    % reference temperature (K)

% set initial condition so zi = 840 m, zb = 610 m, qt = qt0;
zi0 = 840;
zb0 = 610;

% change variables from [zi, zb, qt] to [zi, h, qt];
y0 = MLM_change_vars([zi0;zb0;qtm],'zbq','zhq',SLP,SST,Pref,Tref,do_full,do_simple);

% use alt_y0 for initial conditions if alt_y0 is not empty
if length(alt_y0)==3
   y0 = alt_y0;
end

% Horizontal advection block -- choose hadv(T) and hadv(q) to keep fixed RH
% at surface
if do_simple
   qsf = @(P,T) qs_simple(P,T);
else
   qsf = @(P,T) qs(P,T);
end
   
if do_hadv
   qs_ref = qsf(Pref,Tref);
   rh     = qtm/qs_ref; % BL mean relative humidity ?
   
   hadvq  = rh*dqsdT(Pref,Tref)*hadvT;
   hadvh  = Cp*hadvT + L*hadvq;
else
   hadvq  = 0;
   hadvh  = 0;
end
hadvzi = 0; % for now

%% Initial profile
[dydt,zb,T0,Tl,zref,qs_surf,hs_surf,dqdzu,dqdzs,rref,ep,mu,B,sv0,svl,...
   ql,sv,dqi,dhi,dsvi,we, W0, Wi, E0, Ei, wstar, BIR, rho,...
   wb_bar, wsv_bar, wh_bar, wq_bar,zgrid,Frad,wsed,driz,dEdz,dWdz,...
   lwp,qv,P, LHF, SHF, p, rhod, ib, iFT,hz,qs0,I0,I1,alpha] = BW97_MLM_rhs_v2(0,y0,qt_trop,hp,D,SST,SLP,CTh,CTm,zmax,Nd,...
   do_sed,do_full,a1,a2,ased,LES_tuned,Junya,Pref,Tref,do_driz,...
   do_simple,rad_opts,hadvh,hadvq,hadvzi,do_fixed_LHF,do_fixed_SHF,...
   0,we_funhandle);

do_init_prof = 1;

if do_init_prof
zrel = zgrid./zi;
subplot(2,2,1);
plot(driz.*days_to_sec,zrel,col,-P.*days_to_sec,zrel,[col,'--']); hold on
ylim([0 1.2])

subplot(2,2,2)
plot(wb_bar,zrel,col); hold on
ylim([0 1.2])

subplot(2,2,3)
plot((ql+qv)./gkg,zrel,col,qv./gkg,zrel,[col,'--']); hold on
ylim([0 1.2])

subplot(2,2,4)
plot(hz,zrel,col); hold on
ylim([0 1.2])
end


%% Time evolve the model
rhs = @(t,y) BW97_MLM_rhs_v2(t,y,qt_trop,hp,D,SST,SLP,CTh,CTm,...
   zmax,Nd,do_sed,do_full,a1,a2,ased,LES_tuned,Junya,Pref,Tref,...
   do_driz,do_simple,rad_opts,hadvh,hadvq,hadvzi,do_fixed_LHF,do_fixed_SHF,...
   0,we_funhandle);
options = odeset('RelTol',1e-4,'AbsTol',1e-6); % 'NonNegative',[1 1 1],

[t,y] = ode15s(rhs,tspan,y0,options);

% Initialize time-series output variables
nt = length(t);
SVL{nt} = NaN(size(zgrid));
QC{nt} =  SVL{nt};
SV{nt} =  SVL{nt};
WB{nt} =  SVL{nt};
WSV_BAR{nt} =  SVL{nt};
WH_BAR{nt} =  SVL{nt};
WQ_BAR{nt} =  SVL{nt};
zg{nt} =  SVL{nt};
FR{nt} =  SVL{nt};
WSED{nt} =  SVL{nt};
DRIZ{nt} =  SVL{nt};
HZ{nt} =  SVL{nt};
QV{nt} =  SVL{nt};
Fsed{nt} =  SVL{nt};
Fp{nt} =  SVL{nt};
P_cb =  NaN(size(t));
P_surf =  NaN(size(t));
Fr0 =  NaN(size(t));
Fr_plus =  NaN(size(t));
dzdt   =  NaN(size(t));
dhdt   =  NaN(size(t));
dqdt   =  NaN(size(t));
qli    =  NaN(size(t));
qlb    =  NaN(size(t));
lwp_cb =  NaN(size(t));
lwp_zb =  NaN(size(t));
A =  NaN(size(t));

%% Post-process time-series vars from model output
zi = y(:,1); h = y(:,2); qt = y(:,3);
rho_w = 1000;
for jj=1:length(t)
   tt = t(jj); yy = y(jj,:);
   [dydt,zb(jj,1),T0(jj,1),Tl(jj,1),zref(jj,1),...
      qs_surf(jj,1),hs_surf(jj,1),dqdzu,dqdzs,rref(jj,1),ep,mu,B,sv0(jj,1),SVL{jj},...
   QC{jj},SV{jj},dqi(jj,1),dhi(jj,1),dsvi(jj,1),we(jj,1), W0(jj,1), Wi(jj,1),...
   E0(jj,1), Ei(jj,1), wstar(jj,1), BIR(jj,1), rho,...
   WB{jj}, WSV_BAR{jj}, WH_BAR{jj}, WQ_BAR{jj},zg{jj},FR{jj},WSED{jj},DRIZ{jj},...
   dEdz(jj,1),dWdz(jj,1),lwp_zb(jj,1),QV{jj},P, LHF(jj,1), SHF(jj,1), p, rhod, ib, iFT,...
   HZ{jj},qs0(jj,1),I0(jj,1),I1(jj,1),alpha(jj,1),lwp,A(jj,1)] ...
      = rhs(tt,yy);

   Fsed{jj} = WSED{jj}.*QC{jj}.*rref(jj)./rho_w;     % in (m/s) 
   Fp{jj} = -DRIZ{jj}.*rref(jj)./rho_w+Fsed{jj};     % precip flux (m/s)
   P_cb(jj,1) = interp1(zg{jj},Fp{jj},zb(jj));       % cloud base precip
   P_surf(jj,1) = Fp{jj}(1);                         % surface precip
   Fr0(jj,1) = FR{jj}(1);                            % surface radiation
   Fr_plus(jj,1) = interp1(zg{jj},FR{jj},zi(jj)+50); % Fr(z_i + 50 m)
   dzdt(jj,1) = dydt(1);
   dhdt(jj,1) = dydt(2);
   dqdt(jj,1) = dydt(3);
   qli(jj,1)  = interp1(zg{jj},QC{jj},zi(jj));       % ql at zi-
   qlb(jj,1)  = interp1(zg{jj},QC{jj},zb(jj));       % ql at zb (= 0?)
   lwp_cb(jj,1) = lwp(1);                            % lwp (kg/m3)
end

LHF = (W0 - P_surf).*rref.*L; % (m/s)*rref*L = kg /(m2 s).*(J/kg) = W/m2
SHF = (E0.*rref - Fr0)-LHF;   % sensible heat flux at surface (W/m2)
dFR_BL = Fr_plus - Fr0;       % radiative flux divergence across BL

t = t./days_to_sec; % convert seconds to days

% Plot setup
tmax = t(end);
iNan = find(BIR>0.3,1,'first');
% if ~isempty(iNan)
%    t(iNan:end) = NaN; %don't plot these values
% end

% Steady state using fsolve
[y_fsolve, fval, exitflag, output_fsolve,jac_fsolve] = fsolve(@(z) rhs(t(end)*86400,z),y(end,:));

%% Display summary info
disp(['Nd = ',num2str(Nd)])
disp(['zi = ',num2str(zi(end))])
disp(['zb = ',num2str(zb(end))])
disp(['lwp = ',num2str(lwp_cb(end).*g_per_m2)])
disp(['we = ',num2str(we(end))])
disp(['w* = ',num2str(wstar(end))])
disp(['BIR = ',num2str(BIR(end))])
disp(['P(zb) = ',num2str(P_cb(end).*mm_per_day)])
disp(['P(surf) = ',num2str(P_surf(end).*mm_per_day)])
disp(' ')

%% Make figure 1
figure(11)
subplot(3,1,1)
plot(t,zi,'Color',col); hold on
xlabel('time (days)')
ylabel('z_i (m)')

subplot(3,1,2)
plot(t,h/Cp,'Color',col); hold on
xlabel('time (days)')
ylabel('h/Cp (K)')

subplot(3,1,3)
plot(t,qt*1000,'Color',col); hold on
xlabel('time (days)')
ylabel('q_t (g/kg)')

%% figure 3
figure(13)
subplot(2,3,1)
plot(t,zi,'Color',col); hold on
plot(t,zb,'Color',col,'Linestyle','--');
title('(a) z_i (-) and z_b (--)')
xlabel('t (days)')
ylabel('height (m)')
xlim([0, tmax])

subplot(2,3,2)
plot(t,BIR,'Color',col); hold on
title('(b) Buoyancy Integral Ratio')
xlabel('t (days)')
ylabel('BIR')
xlim([0, tmax])

subplot(2,3,3)
plot(t,lwp_cb.*g_per_m2,'Color',col); hold on
title('(c) Liquid water path')
xlabel('t (days)')
ylabel('LWP (g m^{-2})')
xlim([0, tmax])

subplot(2,3,4)
plot(t,P_cb.*mm_per_day,'Color',col); hold on
plot(t,P_surf.*mm_per_day,'Color',col,'Linestyle','--');
title('(d) cloud base (-) and surf (--) precip.')
xlabel('t (days)')
ylabel('mm day^{-1}')
xlim([0, tmax])

subplot(2,3,5)
plot(t,we.*mm_per_sec,'Color',col); hold on
title('(e) Entrainment rate')
xlabel('t (days)')
ylabel('w_e (mm s^{-1})')
xlim([0, tmax])

subplot(2,3,6)
plot(t,wstar,'Color',col); hold on
title('(f) Convective velocity')
xlabel('t (days)')
ylabel('w^* (m s^{-1})')
xlim([0, tmax])

%% Figure 2            
figure(12)

tvals = (t>t(end)-0.25); ivals = find(tvals);
zbbar = mean(zb(tvals));
zibar = mean(zi(tvals));
zrel  = linspace(0,1.5,151); 
zbrel = zbbar./zibar;

qbar  = mean(qt(tvals))./gkg; % q_t (g/kg)
slbar = mean(h(tvals)-L.*qt(tvals));
qlbbar = mean(qlb(tvals))./gkg;
qlibar = mean(qli(tvals))./gkg;
qt_trop_bar = mean(qt(tvals)+dqi(tvals));
if nargin(hp)==2
 HP = mean(hp(zi(tvals),t(tvals)*86400));
else
 HP = mean(hp(zi(tvals)));
end


xq = [qbar qbar qt_trop_bar./gkg qt_trop_bar./gkg];
subplot(2,3,1)
MLM_profile(zbbar./zibar,1,1.5,xq,'Color',col); hold on
xlabel('q_t (g kg^{-1})')
ylabel('z/z_i')
title('(a) Total water mixing ratio')

xsl = [slbar, slbar,HP-L.*qt_trop_bar, HP-L.*qt_trop_bar]./Cp;
subplot(2,3,2)
MLM_profile(zbbar./zibar,1,1.5,xsl,'Color',col); hold on
xlabel('s_l Cp^{-1} (K)')
ylabel('z/z_i')
title('(b) Liquid static energy')

xql = [0, qlbbar,qlibar, 0, 0];
subplot(2,3,3)
MLM_profile(zbbar./zibar,1,1.5,xql,'Color',col); hold on
xlabel('q_l (g kg^{-1})')
ylabel('z/z_i')
title('(c) Liquid water content')

xdriz = zeros(size(zrel));
xprec = zeros(size(zrel));
xwb   = zeros(size(zrel));

for kk=1:length(ivals)
    ix = ivals(kk);
    xdriz = xdriz + interp1(zg{ix}./zi(ix),-DRIZ{ix}.*rref(ix)./rho_w,zrel);
    xprec = xprec + interp1(zg{ix}./zi(ix),Fp{ix},zrel);
    xwb   = xwb   + interp1(zg{ix}./zi(ix),WB{ix},zrel);
end
xdriz = xdriz ./length(ivals);
xprec = xprec ./length(ivals);
xwb   = xwb ./length(ivals);

subplot(2,3,4)
zplot(zbrel,1.0,xdriz.*mm_per_day,zrel,'Color',col); hold on
zplot(zbrel,1.0,xprec.*mm_per_day,zrel,'Color',col,'Linestyle','--'); hold on
xlabel('precipitation flux (mm day^{-1})')
ylabel('z/z_i')
title('(d) Downward drizzle (-) and precip (--) flux')

subplot(2,3,5)
zplot(zbrel,1.0,xwb,zrel,'Color',col); hold on
xlabel('<w b> (m^2 s^{-3})')
ylabel('z/z_i')
title('(e) Buoyancy flux')
