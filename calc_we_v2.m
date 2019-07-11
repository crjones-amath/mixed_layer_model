function [we, W0, Wi, E0, Ei, wstar, BIR,...
   wb_bar, wsv_bar, wh_bar, wq_bar,Frad,wsed,driz,lwp_zb,P,LHF,SHF,I0,I1,alpha,lwp,A] = ...
   calc_we_v2(zb,zi,zgrid, ql, rho, rref, Nd, sv0, dsvi, hm, qt, dhi, dqi,...
   hsat, qsat, mu, B, ep, CTh, CTm, D, do_sed, ib, iFT, ...
   a1, a2, ased,p,LES_tuned,do_full,Junya,do_driz, do_simple, rad_opts,varargin)
   % calculate the entrainment velocity, radiation, and precipitation flux
   
   %---------------------------------------
   % Input variables:
   %---------------------------------------
   % zb : cloud base
   % zi : inversion
   % zgrid: vertical grid (assumes zi, zb occur at grid points)
   % ql  : ql profile, same length as zgrid
   % rho : air density.  Can either be constant reference density, or a
   %       vector of length zgrid
   % rref : reference density
   % Nd : Droplet concentration
   % sv0: Reference sv
   % dsvi: sv(zi+) - sv(zi-) virtual static energy buoyancy inversion jump
   % hm: h in mixed layer
   % qt: qt in mixed layer
   % dhi: Inversion jump in h;  dhi = hp(zi)  - hm;
   % dqi: Inversion jump of qt; dqi = qt_trop - qt;
   % hsat: h at SST, q = saturation
   % qsat: qt at SST, q = saturated
   % mu: mu from BW97
   % B: beta from BW97
   % ep: epsilon from BW97
   % CTh, CTm: either C*V (surface transfer coeff * surface velocity) for
   %     calculating SHF, LHF, or -- if do_fixed_{LHF,SHF}, CTh = SHF, CTm =
   %      LHF
   % D: Divergence
   % do_sed: 1 do sedimentation; 0 no sedimentation
   % ib = index of zgrid corresponding to zb.  zgrid(ib) = zb;
   % iFT = lowest index of zgrid corresponding to free trop.  
   %       zgrid(iFT-1)= zi
   % a1,a2,ased: TN87 + BBU(2010) entrainment coeffiencts
   % p: pressure.  Either reference pressure or variable
   % LES_tuned: 
   % do_full: 
   % Junya: 
   % do_driz: 
   % do_simple:
   % rad_opts: Structure for passing additional arguments in
   % varargin variables:
   %   do_fixed_LHF = 0
   %   do_fixed_SHF = 0
   %   debug = 0
   %---------------------------------------
   % Output variables:
   %---------------------------------------
   % we: entrainment rate 
   % W0,E0: W(z=0) or E(z=0) from BW97
   % Wi,Ei: W(z=zi-), E(z=zi-) from WB97
   % wstar: Convective velocity scale
   % BIR: Buoyancy integral ratio
   % wb_bar: <w'b'>(z) on zgrid
   % wsv_bar: <w'sv'>(z) on zgrid
   % wh_bar: <w'h'>(z) on zgrid
   % wq_bar: <w'qt'>(z) on zgrid
   % Frad: Radiative flux on zgrid
   % wsed: cloud top sedimentation velocity
   % driz: Drizzle flux on zgrid
   % lwp_zb: column-integrated liquid water path
   % P: Total precipitation flux on zgrid
   % LHF,SHF: Latent and sensible heat flux at surface
   % I0,I1: I_0 and I_1 from BW97 appendix
   % alpha: alpha from BW97 appendix
   % lwp: lwp integrated downward, as a function of z on zgrid
   % A: Entrainment efficiency
   global g L gam
   
   % varargin defaults
   do_fixed_LHF = 0;
   do_fixed_SHF = 0;
   debug = 0;
   
   numvarargs = length(varargin);
   optargs = {do_fixed_LHF,do_fixed_SHF,debug};
   optargs(1:numvarargs) = varargin;
   [do_fixed_LHF, do_fixed_SHF,debug] = optargs{:};
   
   do_lin = ~do_full;
   qc = ql;
   qcmax = ql(iFT-1); % qc at inversion
   pri = p(iFT-1); % needed?!

lwp = LWP(zgrid,qc,rho); % in kg/m2; rho = rho_dry or rref
lwp_zb = lwp(1); % LWP = constant below cloud base, so choose surf value
% lwp_zb = interp1(zgrid,lwp,zb-10); % LWP = constant below cloud base

%-------------------------------------------
% Call Radiation Code
%-------------------------------------------

% Simple radiation:  Get parameters from rad_opts, if available
Frad = zeros(size(zgrid)); % fill temp;
if ~isempty(rad_opts)
    if rad_opts.dosw
       Frad = Frad + sw_simple(zgrid,qc,rho,zi,rad_opts.F0sw, rad_opts.ksw, rad_opts.b0, rad_opts.m0, rad_opts.m1);
    end
    if rad_opts.dolw
       Frad = Frad + rad_simple(zgrid,qc,rho,zi,rad_opts.coef3,rad_opts.F0,rad_opts.F1,rad_opts.k);
    end
    if rad_opts.cldtop_dFR > 0 % concentrated at cloud top, no structure below
        Frad = 0.*Frad;
        Frad(iFT:end) = rad_opts.cldtop_dFR;
    end
else
   Frad = rad_simple(zgrid,qc,rho,zi); % simplest scheme being used ATM
end
Fr0  = Frad(1); %W/m2
Frp = interp1(zgrid(iFT:end),Frad(iFT:end),zi+50); % Fr+ = Fr 50 m above inversion - should i revisit this?


%-------------------------------------------
% Call Sedimentation and Precipitation Code
%-------------------------------------------
[wsed, P] = sedi_flux(qc,Nd,rho); % use simplest scheme for now
if ~do_sed
   wsed = 0.*wsed;
   P = 0.*P;
end

if Junya
   driz = do_driz.*drizzle_flux(zgrid,zb,zi,lwp_zb,Nd,LES_tuned,32.0,Junya); % Junya's version
   Fprec = (driz-P); % convert wsed to moisture precip flux in m/s
else
   driz = do_driz.*drizzle_flux(zgrid,zb,zi,lwp_zb,Nd,LES_tuned); % mm/s = (kg / m2 s)
   Fprec = driz./rref-P./rref; % convert wsed to moisture precip flux in m/s
end

wsedi = wsed(iFT-1); % wsed at zi-
Fp = @(z) interp1(zgrid(1:iFT-1),Fprec(1:iFT-1),z); % Precipitation flux
Fr = @(z) interp1(zgrid(1:iFT-1),Frad(1:iFT-1),z); % radiation flux

%---------------------
% inversion jumps
%---------------------
if do_fixed_LHF
   LHF = CTm; % overload CTm = LHF_in, since both will not be used at once
   W0  = LHF./(rref.*L) + Fprec(1);
else
   W0 = CTm.*(qsat - qt) + Fprec(1);
   LHF = (W0 - Fprec(1)).*rref.*L;
end

if do_fixed_SHF
   SHF = CTh;
   E0  = (SHF + LHF + Fr0)./rref;
else
   E0  = CTh.*(hsat - hm) + Fr0./rref; % Eq. (6)
   SHF = (E0.*rref - Fr0)-LHF;
end

zSC = zgrid(1:ib);     % subcloud grid
zC  = zgrid(ib:iFT-1); % cloud grid.  zgrid(iFT-1) = zi. zgrid(ib) = zb.

% calculate c0, c1 coefficients in the subcloud (SC) and cloud (C) layers.
[c0_SC, c1_SC] = c_calc(zSC,zb,zi,dhi,dqi,E0,W0,Frad(1:ib),Frp,Fprec(1:ib),sv0,B,ep,mu,rref,1);
[c0_C, c1_C]   = c_calc(zC,zb,zi,dhi,dqi,E0,W0,Frad(ib:iFT-1),Frp,Fprec(ib:iFT-1),sv0,B,ep,mu,rref,0);

% integrate to get I0, I1.
if zb>=zi
   I0c = 0;
   I1c = 0;
else
   I0c = trapz(zC,c0_C);
   I1c = trapz(zC,c1_C);
end

if zb<=0
   I0sc = 0;
   I1sc = 0;
else
   I0sc = trapz(zSC,c0_SC);
   I1sc = trapz(zSC,c1_SC);
end
I0 = I0sc + I0c; % (A3)
I1 = I1sc + I1c; % (A3)

alpha = -1;

% instead, use simplified we specifications.
if isfield(rad_opts,'wefixed') % fixed entrainment
   we = rad_opts.wefixed;
   A  = NaN;
elseif isfield(rad_opts,'we_propto_zc') % we = c*(zi-zb)
   we = (rad_opts.we_propto_zc).*(zi-zb);
else

if debug
   subplot(2,1,1)
   plot(c0_SC,zSC,'b',c0_C,zC,'b');
   title(['I0 = ',num2str(I0)]);
   xlabel('c_0')
   ylabel('z')
   
   subplot(2,1,2)
   plot(c1_SC,zSC,'b',c1_C,zC,'b');
   title(['I1 = ',num2str(I1)]);
   xlabel('c_1')
   ylabel('z')   
end

db     = g./sv0.*dsvi;       % delta b, buoyancy jump
dbs    = B.*(dhi./hm - dqi); % delta b_s

if do_lin
   xis = qcmax./(-dqi + (gam./(1 + gam))*(dhi/L));
else
%   xis = qcmax./(-dqi + (gam./(1 + gam))*(dhi/L));
   xis = xi_sat(hm+dhi,hm,qt+dqi,qcmax,qt,zi,pri,do_simple);
end

if Junya
   db = dsvi;
   Cp_delta2 = B.*dhi - ep.*L.*dqi;
   chistar = xis;
   dsvm = 2*(   .5*   chistar *(chistar*Cp_delta2)+ (.5)*(1-chistar)*(dsvi+chistar*Cp_delta2));
   dbs = dsvm;
   xis = 1; % in Junya's version, chistar is absorbed into 1-dsvm/dsvi;
end

if (~do_sed || ased == 0)
    % entrainment efficiency can be calculated explicitly
   A = entr_eff(a1,a2,0,xis,dbs,db);
else
   % Iteration required to find entrainment efficiency
   wstar0 = 1.0; % initial guess
   [wstar, fval, flag] = fzero(@wstar_sedi_closure,wstar0,...
         optimset('TolX',1e-6,'TolFun',1e-12));
   A = entr_eff(a1,a2,ased,xis,dbs,db,wsedi,wstar);
end

alpha = 2.5.*A.*sv0./(g.*zi.*dsvi); % (A3)
we = alpha.*I0 ./ (1 - alpha.*I1); % ENTRAINMENT RATE!!

end

if we<0
   warning('calc_we:neg_we','Entrainment rate negative!')
elseif isnan(we)
   warning('calc_we:NaN_we','Entrainment rate is NaN!');
end
      
Ei = -we.*dhi + Frp./rref;
Wi = -we.*dqi;

% calculate diagnostic MLM buoyancy flux
[wb_bar, wsv_bar, wh_bar, wq_bar, wstar, BIR] = MLM_buoyflux(zgrid,...
    zi,zb,W0,Wi,E0,Ei,mu,L,B,ep,g,rref,sv0,Frad,Fprec);
 
   % for full entrainment closure including sedimentation, A needs to
   % be determined by solving a nonlinear equation
   % wstar^3 = 2.5*[I0 + A(wstar)./(g zi).*sv0./dsvi.*I1];
   % Note: nested function (think carefully about variable scope here)
   function res = wstar_sedi_closure(win)
      % returns residiual (iterate to find consistent A)
      A    = entr_eff(a1,a2,ased,xis,dbs,db,wsedi,win);
      kap  = ( A.*sv0./(g.*zi.*dsvi) ).*2.5;
      we1  = kap.*I0./(1-kap.*I1);
      wout = (2.5*we1./kap).^(1/3);
      res  = wout-win;
   end

end

function [c0, c1] = c_calc(z,zb,zi,dh,dq,E0,W0,Fr,Frp,Fp,sv0,B,ep,mu,rref,isSC)
global g L
if isSC
   zSC_test = (z<=zb);
   zC_test  = (z>zb & z<=zi);
else
   zSC_test = (z<zb);
   zC_test  = (z>=zb & z<=zi);
end
ch =  g./sv0.*(     1.*zSC_test +     B.*zC_test );
cq = -g./sv0.*( mu.*L.*zSC_test + ep.*L.*zC_test );

zs = z./zi; % zeta, scaled z

c0 = ch.*( (1-zs).*E0 - Fr./rref + zs.*Frp./rref) + ...
     cq.*( (1-zs).*W0 - Fp);
c1 = -zs.*(ch.*dh + cq.*dq);
end