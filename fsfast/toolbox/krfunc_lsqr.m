function z = krfunc_lsqr(w,rimgsize,kimgsize,coilprof,transpstring)
% z = krfunc_lsqr(w,rimgsize,kimgsize,coilprof,transpstring)


%
% krfunc_lsqr.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
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

z = [];
if(nargin < 4 | nargin > 5)
  fprintf('z = krfunc_lsqr(w,rimgsize,kimgsize,coilprof,transpstring)\n');
  error('');
  return;
end

if(nargin == 4)
  % No transpose, then run the forward model
  % where w = reconned image, and z = kspace image
  %fprintf('kimg:\n');
  z = kimgfunc(w,rimgsize,kimgsize,coilprof);
else
  % Transpose, then run the backprojection model
  % where w = kspace image, and z = reconned image
  %fprintf('rimg:\n');
  z = rimgfunc(w,rimgsize,kimgsize,coilprof);
end

return;
  
  





