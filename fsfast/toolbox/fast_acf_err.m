function acferr = fast_acf_err(racf,R,yacf)
% acferr = fast_acf_err(racf,R,yacf)


%
% fast_acf_err.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:03 $
%    $Revision: 1.3 $
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

if(nargin ~= 3)
  fprintf('acferr = fast_acf_err(racf,R,yacf)\n');
  return;
end

nf = length(racf);
Vy = toeplitz(yacf);

% Compute expected racf %
for l = 1:nf
  Dl = diag(ones(nf-l+1,1),l-1);  
  racfexp(l) = trace(R*Dl*R*Vy);
end
racfexp = racfexp/racfexp(1);




