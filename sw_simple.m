function F = sw_simple(z,qc,rho,zi,varargin)
% Idealized radiation scheme for SW component.
% Includes exponential decay within the cloud, linear 
% decreases above and below the cloud (assumption: constant absorbtion)

% z   = vertical grid
% qc  = cloud liquid water at each z grid point (size = size(z))
% rho = density of dry air
% optional inputs:
%   varargin(1) = F0
%   varargin(2) = k
%   varargin(3) = b0
%   varargin(4) = m0
%   varargin(5) = m1

% defaults: (start by eyeballing these badboys
F0 = 17; 
k  = 85;
b0 = 164;
m0 = (168-164)/572;
m1 = (7/750) - m0;

% load varargin values
numvarargs = length(varargin);
if numvarargs>5
   error('rad_simple:TooManyArgs','Expect at most 5 optional arguments')
end

% set defaults for optional inputs
optargs = {F0, k, b0, m0, m1};
optargs(1:numvarargs) = varargin;
[F0, k, b0, m0, m1] = optargs{:};

%--------------------------
% begin calculation
%--------------------------
Qc = LWP(z,qc,rho);
Fcloud = F0.*exp(-k.*Qc);

% below (zb, zi?)
Fsurf = b0 + m0.*z; % linear from surface upward
Ftrop = m1.*(z-zi); % adjust the slope above zi if needed

F = (-1)*( Fcloud + Ftrop.*(z>=zi) + Fsurf );
