function zplot(zb, zi, varargin)
ish = ishold;
plot(varargin{:}); hold on
% plot(xlim,[zb zb],'k:',xlim,[zi zi],'k:');
if ~ish
   hold off
end
end
