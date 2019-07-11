function [zgrid, zSC, zC, zFT, zg2] = MLM_grid(zb,zi,zrange,npts)
% Create a grid suitable for use with this MLM, which is 
% essentially a grid with a subcloud (SC) layer, cloud (C) layer,
% and free troposphere (FT) layer, separated by zb and zi, respectively.  
% zb and zi are gridpoints in zgrid, and the grid is non-uniform
% ----------------------------------------------------------------
% Input:
% zb = cloud base (separates SC from C layer)
% zi = inversion  (separates C from FT layer)
% zrange = range of z values for grid, either [zmin, zmax] or just
%      zmax (in which case zmin = 0).
% npts = points in zgrid.  If npts is scalar, each sublayer has npts
%        If npts is given as a vector, then npts = [nSC, nC, nFT].
% -----------------------------------------------------------------
% Output:
%  zgrid = full grid spanning [zmin, zmax]
%  zSC   = SC grid [zmin, zb]
%  zC    = C  grid [zb, zi]
%  zFT   = FT grid [zi, zmax]
%-----------------------------------------------------------------
if length(zrange) == 1
   zmin = 0; zmax = zrange;
elseif length(zrange) == 2
   zmin = zrange(1); zmax = zrange(2);
else
   error('MLM_grid:zrange','zrange should have 2 or fewer elements')
end

if isempty(npts)
   n_dft = 50;
   nSC = n_dft; nC = n_dft; nFT = n_dft;
elseif length(npts)==1
   nSC = npts; nC = npts; nFT = npts;
elseif length(npts)==3
   nSC = npts(1); nC = npts(2); nFT = npts(3);
else
   error('MLM_grid:npts','npts must be scalar or vector of length 3')
end

if zb>=zi
   warning('MLM_grid:no_cloud','zb>zi, so no cloud (stuff may break)')
   zSC = linspace(zmin,zi,nSC);
elseif zb<=zmin
   warning('MLM_grid:no_subcloud','zb<=zmin, so only fog (stuff may break)')
   zSC = linspace(zmin,zmin,nSC);

else
   zSC   = linspace(zmin, zb, nSC); % subcloud z grid
end

zC    = linspace(zb,zi, nC);     % cloudy   z grid
zFT   = linspace(zi,zmax, nFT);  % free troposphere z grid
zgrid = unique([zSC,zC,zFT]);    % one big happy grid
zg2   = [unique([zSC,zC]),zFT];  % one big happy grid with 2 entries for zi (one above, one below)
end