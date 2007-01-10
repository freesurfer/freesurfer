function [w, Sigma_s, Sigma_r, M] = fast_wopt(y,Xf,Xr)
% [w, Sigma_s, Sigma_r] = fast_wopt(y,Xf,<Xr>)


%
% fast_wopt.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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

w = [];
Sigma_s = 0;
Sigma_r = 0;
if(nargin < 2 | nargin > 3)
  fprintf('[w, Sigma_s, Sigma_r] = fast_wopt(y,Xf,<Xr>)\n');
  return;
end

if(~exist('Xr','var')) Xr = []; end


X = [Xf Xr];
Xf0 = [Xf zeros(size(Xr))];

[ntp nobs] = size(y);
dof = ntp - size(X,2);
if(dof < nobs )
  fprintf('ERROR: dof = %d < number of observations = %d\n',...
	  dof,nobs);
  return;
end

beta = (inv(X'*X)*X')*y;
r = y - X*beta;
s = Xf0*beta;

Sigma_r = r'*r;
Sigma_s = s'*s;

M = Sigma_s*inv(Sigma_r);

[u s v] = svd(M);

w = u(:,1);

return;





  
  
  
