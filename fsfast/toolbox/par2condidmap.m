function condidmap = par2condidmap(par,nullid)
%
% condidmap = par2condidmap(par,nullid)
%
% Extracts and sorts all unique condition numbers in paradigm.
%
%


%
% par2condidmap.m
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

if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: condidmap = par2condidmap(par,<nullid>)'
  qoe(msg); error(msg);
end

condid = reshape1d(par(:,2,:,:));
condidmap = unique(condid);

if(nargin == 2)
  n = find(condidmap == nullid);
  if(isempty(n))
    condidmap = [nullid; reshape1d(condidmap)];
  else
    tmp = condidmap(1);
    condidmap(1) = nullid;
    condidmap(n) = tmp;
  end
end


return;
