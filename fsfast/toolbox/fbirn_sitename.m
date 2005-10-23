function sitename = fbirn_sitename(siteno)
% sitename = fbirn_sitename(siteno)
%
% Returns the site name given the site number.
% 
% $Id: fbirn_sitename.m,v 1.1 2005/10/23 19:39:42 greve Exp $

if(nargin ~= 1)
  fprintf('sitename = fbirn_sitename(siteno)\n');
  return;
end

sitename = [];
switch(siteno)
 case 1, sitename = 'dunc'; 
 case 2, sitename = 'iowa'; 
 case 3, sitename = 'nm'; 
 case 4, sitename = 'uci'; 
 case 5, sitename = 'ucsd'; 
 case 6, sitename = 'bwh';
 case 7, sitename = 'mgh'; 
 case 8, sitename = 'min'; 
 case 9, sitename = 'stan'; 
 case 10, sitename = 'dunc4t';
end

return;  
  






