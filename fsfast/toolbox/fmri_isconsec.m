function r = fmri_isconsec(seq,nconsec)
%
% r = fmri_isconsec(seq,nconsec)
%
% Returns 1 if any stimulus in seq is consecutively
% presented more than nconsec times.  Returns 0 if not.
%
%


%
% fmri_isconsec.m
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

if(nargin ~= 2)
  msg = 'Usage: r = fmri_isconsec(seq,nconsec)';
  qoe(msg);error(msg);
end

r = 1;

if(nconsec == 1) dtest = 1;
else             dtest = 0;
end

for n = min(seq):max(seq),
  ind = find(seq==n);
  for m = 1:nconsec,
    d = diff(ind,1);
    if(length(d)==0) break; end
    ind = find(d==1);
  end
  if(~isempty(d))
    if(isempty(find(d==1))) return; end
  end
end

%  ind = find(seq==n);
%  d = diff(ind,nconsec);
%  if(~isempty(find(d==dtest))) return; end


r = 0;
return;
