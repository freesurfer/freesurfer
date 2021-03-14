function sitename = fbirn_sitename(siteno)
% sitename = fbirn_sitename(siteno)
%
% Returns the site name given the site number.
% 
%


%
% fbirn_sitename.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

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
  






