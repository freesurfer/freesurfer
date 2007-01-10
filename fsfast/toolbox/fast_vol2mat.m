function mat = fast_vol2mat(vol,sliceflag)
% mat = fast_vol2mat(vol,<sliceflag>)
%
% Reshapes a volume (size [nr nc ns nf]) into a matrix
% of size nf by (nr*nc*ns). vol can also be an mri struct.
% See MRIread.m.
%
% [nr nc ns nf] = size(vol).
%
% If sliceflag=1, then [nr nc nf] = size(vol), ie, the volume
% is really an image, and the 3rd dim is treated as the frame
% dimension instead of the 3rd spatial dimension.
%
% See also: fast_mat2vol.
%
%


%
% fast_vol2mat.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
%    $Revision: 1.4 $
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

mat = [];
if(nargin < 1 | nargin > 2)
  fprintf('mat = fast_vol2mat(vol,<sliceflag>)\n');
  return;
end

if(~exist('sliceflag','var')) sliceflag = []; end
if(isempty(sliceflag)) sliceflag = 0; end

if(isfield(vol,'vol'))  vol = vol.vol; end

if(length(size(vol)) > 3 & sliceflag) 
  fprintf('ERROR: sliceflag set but volume is 4D\n');
  return;
end

[nr nc ns nf] = size(vol);
nv = nr*nc*ns;
if(sliceflag)
  nf = ns;
  nv = nr*nc;
end

mat = transpose(reshape(vol,[nv nf]));


return;
