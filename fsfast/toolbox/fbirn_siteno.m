function siteno = fbirn_siteno(sitename)
% siteno = fbirn_siteno(sitename)
%
% Returns the site number id given the site name
% 
% $Id: fbirn_siteno.m,v 1.2 2005/10/23 21:05:09 greve Exp $

if(nargin ~= 1)
  fprintf('siteno = fbirn_siteno(sitename)\n');
  return;
end

siteno = [];

sitelist = fbirn_sitelist;
nsites = size(sitelist,1);
for nthsite = 1:nsites
  if(strcmp(sitename,deblank(sitelist(nthsite,:))))
    siteno = nthsite;
    return;
  end
end


return;  
  






