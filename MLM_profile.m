function MLM_profile(zb,zi,zmax,x,varargin)
% plot MLM profile of x, where x is piecewise linear b/w 0,zb,zi,zmax
% zi   = inversion
% zb   = cloud base
% zmax = max plot
% x    = value to be plotted, 
%   w/ x = [x(0) x(zi-) x(zi+) x(zmax)]
%   or x = [x(0) x(zb) x(zi-) x(zi+) x(zmax)]
%   or x = [x(0) x(zb-) x(zb+) x(zi-) x(zi+) x(zmax)]
% varargin = additional plot options
n = length(x);
do_hold = ishold; % test if current plot has hold on
hold on; % turn hold on for the duration of this function; return to do_hold after

if n== 4
   zvals = [0 zi zi zmax];
elseif n==5
   zvals = [0 zb zi zi zmax];
elseif n==6
   zvals = [0 zb zb zi zi zmax];
else
   error('MLM_profile:xvals','length(x) must be 5 or 6')
end
plot(x,zvals,varargin{:});
%plot(xlim,[zb zb],'k:',xlim,[zi zi],'k:');

if ~do_hold
   hold off
end
end