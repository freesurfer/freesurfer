function p = mri_cdf2p(x,xcdf,cdf)
% p = mri_cdf2p(x,xcdf,cdf)
%
% Prob that x <= X. x can be any size, but xcdf and cdf
% must be vectors of the same length.
%
% Note: when creating xcdf and cdf from ComputePDF,
% make sure to subtract xdelta from xcdf before
% using this function.
%
% $Id: mri_cdf2p.m,v 1.1 2008/02/05 00:05:32 greve Exp $

%
% mri_cdf2p.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2008/02/05 00:05:32 $
%    $Revision: 1.1 $
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

p = [];
if(nargin ~= 3)
  fprintf('p = mri_cdf2p(x,xcdf,cdf)\n');
  return;
end

xsize = size(x);
x = transpose(x(:)); % row vector
nx = length(x);

xcdf = xcdf(:); % column vector
nxcdf = length(xcdf);

dx = repmat(x,[nxcdf 1]) - repmat(xcdf,[1 nx]);

[dxmin idxmin] = min(abs(dx));

p = cdf(idxmin);
p = reshape(p,xsize);

return;


