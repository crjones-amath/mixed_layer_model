function vout = MLM_change_vars(vin,in_type,out_type,p0,sst,Pref,Tref,do_full,do_simple)
global L
% transform between sets of thermodynamic variables in MLM
%------------------------------------------------------------
% types:
% 1.  v = [zi, h    , qt] (zhq)
% 2.  v = [zi, zb   , sl] (zbs)
% 3.  v = [zi, zi-zb, sl] (zcs)
% 4.  v = [zi, sl   , qt] (zsq)
% 5.  v = [zi, zb   , qt] (zbq)
% 6.  v = [zi, zb   ,svl] (zbv)
%------------------------------------------------------------
% varargin:
% do_full   = {1 if using full thermo, 0 if lin thermo}
% do_simple = (1 if using qs_simple, 0 if using qs)
%------------------------------------------------------------
% Set defaults for variable input arguments
% numvarargs = length(varargin);
% do_full   = 1;
% do_simple = 0; 
% optargs = {do_full, do_simple};
% optargs(1:numvarargs) = varargin;
% [do_full, do_simple] = optargs{:};

zb_fun = @(h,q)   zb_from_h_qt(h,q,p0,sst,Pref,Tref,do_full, do_simple);
qt_fun = @(sl,zb) qt_from_sl_zb(sl,zb,p0,sst,Pref,Tref,do_full, do_simple);
sl_fun = @(q,zb) sl_from_qt_zb(q,zb,p0,sst,Pref,Tref,do_full, do_simple);
qt_fun2 = @(svl,zb) qt_from_svl_zb(svl,zb,p0,sst,Pref,Tref,do_full, do_simple);

%-----------------------------------------
% (calculate missing thermo variables from vin
%-----------------------------------------
zi = vin(1); %same for all
mu = muL(Tref)./Lv(Tref); % only needed for case{6}

switch in_type
   case {1,'zhq'}
      h  = vin(2); qt = vin(3);
      sl = h - L.*qt;
      zb = zb_fun(h,qt);
      svl = h - mu.*L.*qt;
   case {2,'zbs'}
      zb = vin(2); sl = vin(3);
      qt = qt_fun(sl,zb);
      h  = sl + L.*qt;
      svl = h - mu.*L.*qt;
   case {3,'zcs'}
      zb = zi-vin(2); sl = vin(3);
      qt = qt_fun(sl,zb);
      h  = sl + L.*qt;      
      svl = h - mu.*L.*qt;
   case {4,'zsq'}
      sl = vin(2); qt = vin(3);
      h  = sl + L.*qt;
      zb = zb_fun(h,qt);
      svl = h - mu.*L.*qt;
   case {5,'zbq'}
      zb = vin(2); qt = vin(3);
      sl = sl_fun(qt,zb);
      h  = sl + L.*qt;
      svl = h - mu.*L.*qt;
   case {6,'zbv'}
      zb = vin(2); svl = vin(3);      
      qt = qt_fun2(svl,zb);
      h  = svl + mu.*L.*qt;
      sl = h - L.*qt;
   otherwise
      error('MLM_change_vars:in_type','Unknown in_type')
end

if qt<0
   warning('MLM_change_vars:qt_neg','qt < 0')
end

if zb>zi
   warning('MLM_change_vars:zb_toolarge','zb > zi')
end

%--------------------------------------
% choose correct output format
%--------------------------------------
vout = vin; % to keep the same shape
switch out_type
   case {1,'zhq'}
      vout(2) = h; vout(3) = qt;
   case {2,'zbs'}
      vout(2) = zb; vout(3) = sl;
   case {3,'zcs'}
      vout(2) = zi-zb; vout(3) = sl;
   case {4,'zsq'}
      vout(2) = sl; vout(3) = qt;
   case {5,'zbq'}
      vout(2) = zb; vout(3) = qt;
   case {6,'zbv'}
      vout(2) = zb; vout(3) = svl;
   otherwise
      error('MLM_change_vars:out_type','Unknown out_type')
end

end

function zb = zb_from_h_qt(h,qt,p0,sst,Pref,Tref,do_full, do_simple)
global L Cp Rd delta g Rv

sl = h - L.*qt; % may be useful

if do_simple
   qsf = @(p,T) qs_simple(p,T);
else
   qsf = @(p,T) qs(p,T);
end

if do_full
   pb = @(z) p0.*((sl-g.*z)./sl).^(Cp./(Rd.*(1+delta.*qt)  ) );
   Tb = @(z) 1./Cp.*(sl-g.*z);
   
   zb = fzero(@(z) qsf(pb(z),Tb(z)) - qt, 400);
else
   dqdzu = dqsdzu(Pref,Tref);
   if do_simple
      qs_ref  = qsf(Pref,Tref);
      dqsdT   = qs_ref.*L./(Rv.*Tref.^2);
      dqsdp   = -qs_ref./Pref;
      gam     = (L./Cp).*dqsdT;
      H_scale = Rd.*Tref./g;
      dqdzu   = (  (Rd.*Tref./Cp) .* dqsdT + Pref.*dqsdp)./H_scale;

      qs_surf = qsf(p0,sst);
      hs_surf = Cp.*sst + L.*qs_surf;
      
      zb = ((1+gam).*(qs_surf-qt) - gam./L.*(hs_surf-h))./dqdzu;
   else
      qs0 = qsf(p0,sl./Cp); 
      zb = max(0, (qs0 - qt)./dqdzu);
   end
end

end

function qt = qt_from_sl_zb(sl,zb,p0,sst,Pref,Tref,do_full, do_simple)
global L Cp Rd delta g Rv

if do_simple
   qsf = @(p,T) qs_simple(p,T);
else
   qsf = @(p,T) qs(p,T);
end

if do_full
   pb = @(q) p0.*((sl-g.*zb)./sl).^(Cp./(Rd.*(1+delta.*q)  ) );
   Tb = 1./Cp.*(sl-g.*zb);
   
   qt = fzero(@(q) qsf(pb(q),Tb) - q, 8e-3);
else
   dqdzu = dqsdzu(Pref,Tref);
   if do_simple
      qs_ref  = qsf(Pref,Tref);
      dqsdT   = qs_ref.*L./(Rv.*Tref.^2);
      dqsdp   = -qs_ref./Pref;
      gam     = (L./Cp).*dqsdT;
      H_scale = Rd.*Tref./g;
      dqdzu   = (  (Rd.*Tref./Cp) .* dqsdT + Pref.*dqsdp)./H_scale;

      qs_surf = qsf(p0,sst);
      hs_surf = Cp.*sst + L.*qs_surf;
      
      qt = (1+gam).*qs_surf - dqdzu.*zb - gam./L.*(hs_surf-sl);
   else
      qs0 = qsf(p0,sl./Cp); 
      qt = qs0 - dqdzu.*zb;
   end
end

end

function sl = sl_from_qt_zb(qt,zb,p0,sst,Pref,Tref,do_full, do_simple)
global L Cp Rd delta g Rv

if do_simple
   qsf = @(p,T) qs_simple(p,T);
else
   qsf = @(p,T) qs(p,T);
end

if do_full
   
   pb = @(sl) p0.*((sl-g.*zb)./sl).^(Cp./(Rd.*(1+delta.*qt)  ) );
   Tb = @(sl) 1./Cp.*(sl-g.*zb);
   
   sl = fzero(@(s) qsf(pb(s),Tb(s)) - qt, 2.96e5);
else
   dqdzu = dqsdzu(Pref,Tref);
   if do_simple
      qs_ref  = qsf(Pref,Tref);
      dqsdT   = qs_ref.*L./(Rv.*Tref.^2);
      dqsdp   = -qs_ref./Pref;
      gam     = (L./Cp).*dqsdT;
      H_scale = Rd.*Tref./g;
      dqdzu   = (  (Rd.*Tref./Cp) .* dqsdT + Pref.*dqsdp)./H_scale;

      qs_surf = qsf(p0,sst);
      hs_surf = Cp.*sst + L.*qs_surf;
      
      sl = -(((1+gam).*qs_surf - dqdzu.*zb - qt).*L./gam - hs_surf);
%      qt = (1+gam).*qs_surf - dqdzu.*zb - gam./L.*(hs_surf-sl);
   else
      sl = fzero(@(s) qsf(p0,s./Cp) - qt - dqdzu.*zb,Cp.*Tref);
%       qs0 = qsf(p0,sl./Cp); 
%       qt = qs0 - dqdzu.*zb;
   end
end

end

function qt = qt_from_svl_zb(svl,zb,p0,sst,Pref,Tref,do_full, do_simple)
global L Cp Rd delta g Rv

mu = muL(Tref)./Lv(Tref);

if do_simple
   qsf = @(p,T) qs_simple(p,T);
else
   qsf = @(p,T) qs(p,T);
end

   slf = @(q) svl - (1-mu).*L.*q; % svl = h - L*qt + L*qt - muLqt = sl + L(1-mu)qt

if do_full
   pb = @(q) p0.*((slf(q)-g.*zb)./slf(q)).^(Cp./(Rd.*(1+delta.*q)  ) );
   Tb = @(q) 1./Cp.*(slf(q)-g.*zb);
   
   qt = fzero(@(q) qsf(pb(q),Tb(q)) - q, 8e-3);
else
   dqdzu = dqsdzu(Pref,Tref);
   if do_simple
      qs_ref  = qsf(Pref,Tref);
      dqsdT   = qs_ref.*L./(Rv.*Tref.^2);
      dqsdp   = -qs_ref./Pref;
      gam     = (L./Cp).*dqsdT;
      H_scale = Rd.*Tref./g;
      dqdzu   = (  (Rd.*Tref./Cp) .* dqsdT + Pref.*dqsdp)./H_scale;

      qs_surf = qsf(p0,sst);
      hs_surf = Cp.*sst + L.*qs_surf;
      
      qt = ((1+gam).*qs_surf - dqdzu.*zb - gam./L.*(hs_surf-svl))./(1+gam-mu.*gam);
   else
      qs0 = @(q) qsf(p0,slf(q)./Cp);
%      qt = qs0 - dqdzu.*zb;
      qt = fzero(@(q) qs0(q)-dqdzu.*zb - q, 8e-3);
   end
end

end