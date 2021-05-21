function [yacf, M] = fast_yacf_kjw(racf,R,p)
% yacf = fast_yacf_kjw(racf,R,<p>)
% 
% Remove bias from residual autocorrelation function using Keith
% Worsley's fmristats method. The bias is induced by projecting out
% the task components to form the residual.
%
% racf is the autocor function of the residual. This
%  can be nf X nv. See fast_acorr.m
% R is the residual forming matrix
% yacf is an approximation of the autocor function
%  of the original noise.
% M is the correction matrix.
% p - cant remember. Default is number of frames. I think it's the
%  same as running it like fast_yacf_kjw(racf(1:p,:),R)
%
% Notes:
%  1. As the number of frames increases, the computation is 
%     more intense and the correction matrix is more singular.
%     There is some happy medium where the accuracy of the 
%     ACF is best, but it depends upon the design among other
%     things.
%
% See also: fast_acorr, fast_kjw_mtx.
%
% Worsely, 2002, NI 15, 1-15.
% 
%


%
% fast_yacf_kjw.m
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

if(nargin < 2 | nargin > 3)
  fprintf('yacf = fast_yacf_kjw(racf,R,<p>)\n');
  return;
end


[nf nv] = size(racf);

if(~exist('p','var')) p = nf; end

M = fast_kjw_mtx(R,p);

yacf = inv(M)*racf(1:p,:);
yacf = yacf./repmat(yacf(1,:),[p 1]);

return;
