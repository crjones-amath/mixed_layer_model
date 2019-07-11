function [wb_bar, wsv_bar, wh_bar, wq_bar, wstar,BIR] = ...
   MLM_buoyflux(z,zi,zb,W0,Wi,E0,Ei,mu,L,B,ep,g,rho,sv0,Frad,Fprec)

wh_bar = (1-z./zi).*E0 + (z./zi).*Ei - Frad./rho;
wq_bar = (1-z./zi).*W0 + (z./zi).*Wi - Fprec;

% assume no turbulence above inversion
wh_bar(z>zi) = 0;
wq_bar(z>zi) = 0;

wsv_bar = wh_bar - mu.*L.*wq_bar;
wsv_bar(z>zb & z<=zi) = B.*wh_bar(z>zb & z<=zi) - ep.*L.*wq_bar(z>zb & z<=zi);
% wsv_bar(z>zi) = 0; % assume no turbulence above
wb_bar = (g./sv0).*wsv_bar;

wstar = (2.5.*trapz(z,wb_bar)).^(1/3);

%-------------------------------------------
% while we're at it, calculate BIR as well
%-------------------------------------------
zSC = z<=zb;
zSI = z<=zi;
wb_neg = wb_bar;
wb_neg(wb_neg>0 | z>=zb) = 0;
wb_pos = wb_bar - wb_neg;
wb_pos(z>zi) = 0;
% wb_neg(wb_neg>0) = 0;
% wb_pos = wb_bar(zSI);
% wb_pos(1:length(wb_neg)) = wb_pos(1:length(wb_neg))-wb_neg;
% wb_pos = wb_bar(zSC)-wb_neg;

% %---------------------------
% % alternative version
% %---------------------------
% IBF = cumtrapz(z,wb_bar);
% ib = find(z<=zb,1,'last'); % last index below cloud base;
% Imax = max(IBF(z<=zb)); % max value below CB corresponds to last positive contribution
% Ineg = Imax-IBF(ib); % negative contribution = this difference
% BIR = Ineg./(IBF(end)-Ineg);

try
   BIR = -trapz(z(zSC),wb_neg(zSC))./trapz(z(zSI),wb_pos(zSI));
catch ME
   ME.identifier
   ME.message
   BIR = NaN;
end