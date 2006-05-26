function r = tdr_rtraj(nrsamples)
% r = tdr_rtraj(nrsamples)
%
% recon space 'trajectory' (returned as a row vector) 
% 
% Starts at -(nrsamples/2-1) (eg, -63 for 128 or -31 for 64)
% Ends at    +nrsamples/2    (eg, +64 for 128 or +32 for 64)
% Increments by 1 each sample.
% Passes through 0 at nrsamples/2.
%
% See also tdr_uniform_phtraj.m
%
% $Id: tdr_rtraj.m,v 1.1 2006/05/26 23:49:11 greve Exp $
%

r = [];
if(nargin ~= 1)
  fprintf('r = tdr_rtraj(nrsamples)\n');
  return;
end

% Starts at -63, passes thru 0 at nc/2, ends at +64
rstart = -(nrsamples/2-1);
rend   = +nrsamples/2;
r = [rstart:rend];

return;
