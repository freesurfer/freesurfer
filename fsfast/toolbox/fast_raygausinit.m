function params = fast_raygausinit(y)
% params = fast_raygausinit(y)
% params = [alpha Rmu Gmu Gstd];
%
%


%
% fast_raygausinit.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:05 $
%    $Revision: 1.5 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
  fprintf('params = fast_raygausinit(y)\n');
  return;
end

if(0)
%md = mode(y);
ysort = sort(y);
ny = length(y);
md = ysort(round(ny/2));
indray = find(y < md);
indgaus = find(y >= md);
end

ymn = mean(y);
indray  = find(y < ymn);
indgaus = find(y > ymn);

alpha = .25;
Rmu = mean(y(indray));
Gmu = mean(y(indgaus));
Gstd = std(y(indgaus));

params = [alpha Rmu Gmu Gstd];

return;
