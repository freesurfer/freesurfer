function z = fast_p2z(p)
% z = fast_p2z(p)
%
% Converts a significance to a z-score. p must be
% between -1 and 1
%
% $Id: fast_p2z.m,v 1.1 2004/05/14 04:23:04 greve Exp $
%

z = [];
if(nargin ~= 1) 
  fprintf('z = fast_p2z(p)\n');
  return;
end

z = erfinv(1-abs(p)) .* sign(p);

return;










