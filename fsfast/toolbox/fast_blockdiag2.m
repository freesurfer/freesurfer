function M = fast_blockdiag2(varargin)
% M = fast_blockdiag2(M1,M2,...)
%
% Assmbles matrix inputs into a block diagonal matrix.
% 
% See also fast_blockdiag
%
%


%
% fast_blockdiag2.m
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

M = [];
if(nargin == 0)
  fprintf('M = fast_blockdiag2(M1,M2,...)\n');
  return;
end

ncols = 0;
nrows = 0;
for n = 1:nargin
  Mn = varargin{n};
  nrows = nrows + size(Mn,1);
  ncols = ncols + size(Mn,2);
end

M = zeros(nrows,ncols);

r1 = 1;
c1 = 1;
for n = 1:nargin
  Mn = varargin{n};
  [nr nc] = size(Mn);
  r2 = r1 + nr - 1;
  c2 = c1 + nc - 1;
  M(r1:r2,c1:c2) = Mn;
  r1 = r2 + 1;
  c1 = c2 + 1;
end

return;







return;
