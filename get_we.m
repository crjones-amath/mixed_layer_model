function we = get_we(t,yin,rhs,varargin)

nw = 21;
numvarargs = length(varargin);
optargs = {nw};
   optargs(1:numvarargs) = varargin;
   [nw] = optargs{:};

   y = yin(:); % because I was silly
% nw = 21;
dum = '[';
for j=1:nw-1
    dum = [dum,'v1,'];
end
dum = [dum,'we] = rhs(t,y);'];

eval(dum);