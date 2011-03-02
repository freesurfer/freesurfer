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
% $Id: mri_cdf2p.m,v 1.2 2011/03/02 00:04:12 nicks Exp $

%
% mri_cdf2p.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.2 $
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


