function sitename = fbirn_sitename(siteno)
% sitename = fbirn_sitename(siteno)
%
% Returns the site name given the site number.
% 
% $Id: fbirn_sitename.m,v 1.2 2005/10/23 21:05:09 greve Exp $

sitename = [];

if(nargin ~= 1)
  fprintf('sitename = fbirn_sitename(siteno)\n');
  return;
end

sitelist = fbirn_sitelist;
nsites = size(sitelist,1);

if(siteno > nsites)
  fprintf('ERROR: siteno = %d > nsites = %d\n',siteno,nsites);
  return;
end

sitename = deblank(sitelist(siteno,:));

return;  
  






