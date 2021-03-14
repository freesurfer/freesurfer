function siteno = fbirn_siteno(sitename)
% siteno = fbirn_siteno(sitename)
%
% Returns the site number id given the site name
% 
%


%
% fbirn_siteno.m
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
  






