function siteno = fbirn_siteno(sitename)
% siteno = fbirn_siteno(sitename)
%
% Returns the site number id given the site name
% 
% $Id: fbirn_siteno.m,v 1.1 2005/10/23 19:39:42 greve Exp $

if(nargin ~= 1)
  fprintf('siteno = fbirn_siteno(sitename)\n');
  return;
end

siteno = [];

switch(sitename)
 case 'dunc', siteno = 1;
 case 'iowa', siteno = 2;
 case 'nm', siteno = 3;
 case 'uci', siteno = 4;
 case 'ucsd', siteno = 5;
 case 'bwh', siteno = 6;
 case 'mgh', siteno = 7;
 case 'min', siteno = 8;
 case 'stan', siteno = 9;
 case 'dunc4t', siteno = 10;
end

return;  
  






