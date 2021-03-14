function z = fast_p2z(p)
% z = fast_p2z(p)
%
% Converts a significance to a z-score. p is assumed to be
% one-sided (if p is signed, then the z gets the sign of p)
%
%  z = sqrt(2)*erfcinv(2*abs(p)) .* sign(p);
% Check:
%  r = randn(100000,1);
%  phat = length(find(r>z))/length(r)
%  phat should be same as p (approx)

%
% fast_p2z.m
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

z = [];
if(nargin ~= 1) 
  fprintf('z = fast_p2z(p)\n');
  return;
end

z = sqrt(2)*erfcinv(2*abs(p)) .* sign(p);


return;










