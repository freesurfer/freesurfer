function slice = fast_ldslice(volid,sliceno)
% slice = fast_ldslice(volid,sliceno)


%
% fast_ldslice.m
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

slice = [];

if(nargin ~= 2)
  msg = 'USAGE: slice = fast_ldslice(volid,sliceno)'
  qoe(msg) ; error(msg);
end

fmt = fast_getvolformat(voldid);
if(isempty(fmt))
  msg = sprintf('Could not determine format of %s',volid);
  qoe(msg) ; error(msg);
end

switch(fmt)

end



return;
