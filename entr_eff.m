function A = entr_eff(varargin)
% calculate the entrainment efficiency for MLM
% A = a1*(1 + a2*chi*(1-db_s./db)exp(-ased*wsed./wstar)

numvarargs = length(varargin);
% default options
a1   = 0.2; a2   = 25; ased = 9;
xi  = 0; dbs = 0; db = 1; wsed = 0; wstar = 1;

optargs = {a1,a2,ased,xi,dbs,db,wsed,wstar};
optargs(1:numvarargs) = varargin;
[a1, a2, ased, xi, dbs, db, wsed, wstar] = optargs{:};

A = a1.*(1 + a2.*xi.*(1-dbs./db).*exp(-ased.*wsed./wstar));

end