function racfexp = fast_racf_exp(R,yacf)
% racfexp = fast_racf_exp(R,yacf)
%
% Computes the expected value of the acf of the residual
% given the acf of y and the residual error forming matrix R.
%
% Based on Worsley, et al, NeuroImage 15, 1-15, 2002.
%


%
% fast_racf_exp.m
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
  fprintf('racfexp = fast_racf_exp(R,yacf)\n');
  return;
end

nf = length(yacf);
Vy = toeplitz(yacf);
RVy = R*Vy;

for l = 1:nf
  Dl = diag(ones(nf-l+1,1),l-1);  
  racfexp(l) = trace(R*Dl*RVy);
end
racfexp = racfexp/racfexp(1);

return;


