function [xmax, ixmax] = maxmd(x)
%
%  [xmax ixmax] = maxmd(x)
%
% Finds the max of the ND array x and
% the indicies of the max.  xmax is
% a scalar and ixmax is a 1-D vector
% of length equal to the number of dimensions in x.
%   
% Douglas N. Greve
% February 3, 1999


%
% maxmd.m
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

if(nargin ~= 1)
     error('USAGE: [xmax ixmax] = maxmd(x)');
end

xdim = size(x);
ndim = length(xdim);

[xmax nxmax] = max(reshape(x,[prod(xdim) 1]));

% fprintf('nxmax = %d\n',nxmax);

ixmax = zeros(1,ndim);
k = nxmax;

for n = ndim : -1 : 2,
   p = prod(xdim(1:n-1));
   j = ceil(k/p);
   ixmax(n) = j;
   % fprintf('n=%d, k=%d, p=%d, j=%d\n',n,k,p,j);
   k = k - (j-1)*p;
end

ixmax(1) = k;

return;

