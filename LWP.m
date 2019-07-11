function lwp = LWP(z,qc,rho)
% calculate LWP given z-grid, qv, and ql
% LWP = 0 at cloud top, positive below
% LWP = \int_z^\infty \rho(z) qc(z) dz
% UNITS:  kg / m^2;

lwp = cumtrapz(z(end:-1:1),rho(end:-1:1).*qc(end:-1:1)); %integrate from z to infty
   % should be zero above cloud top, constant below cloud base
lwp = -1.*lwp(end:-1:1); % put back in proper order
end