function [rmse, racfest] = fast_rphacfrmse(params,racf,R,taper)
% [rmse, racfest] = fast_rphacfrmse(params,racf,R,<taper>)


%
% fast_rphacfrmse.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

if(nargin < 1 | nargin > 4)
  fprintf('[rmse racfest]= fast_rphacfrmse(params,racf,R,<taper>)\n');
  return;
end

alpha = params(1);
rho   = params(2);
T     = params(3);

nf = length(racf);
tau = [0:nf-1]';

nacf = fast_rphacf(alpha,rho,T,nf);
racfest = fast_cvm2acor(R*toeplitz(nacf)*R,1);

if(nargin ~= 4)
  rmse = mean( (racf-racfest).^2 );
  rmseb = mean( (racf-nacf).^2 );
  rmse = rmse + rmseb/4;
else
  rmse = mean( ((racf-racfest).*taper).^2 );
end

return;
