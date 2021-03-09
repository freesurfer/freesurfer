function vol = fast_mat2vol(mat,szvol,sliceflag)
%
% vol = fast_mat2vol(mat,szvol,<sliceflag>)
% 
% Reshapes a 2d matrix mat of size nframes X nvoxels
% to a volume of size [szvol nframes].
%
% If mat is an mri struct (see MRIread.m) with field volmat, 
% then only one arg is needed, and uses:
%      szvol = mat.volsize
%      mat = mat.volmat
%
% szvol = [nrows ncols nslices]. If szvol only has two elements,
% then nslices=1. If sliceflag is set, then nslices=1, and the
% number of spatial voxels is assumed to be nrows*ncols.
%
% See also fast_vol2mat.
%
%


%
% fast_mat2vol.m
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

vol = [];

if(nargin < 1 | nargin > 3)
  fprintf('vol = fast_mat2vol(mat,<szvol>,<sliceflag>)\n');
  return;
end

if(isfield(mat,'volmat'))
  szvol = mat.volsize;
  mat   = mat.volmat;
else
  if(nargin < 2)  
    fprintf('ERROR: fast_mat2vol: need szvol if mat is not an mri struct\n');
    return;
  end
end

if(length(szvol) < 2)
  fprintf('ERROR: fast_mat2vol: szvol dim must be >= 2\n');
  return;
end

if(~exist('sliceflag','var')) sliceflag = []; end
if(isempty(sliceflag)) sliceflag = 0; end

if(sliceflag)  szvol = szvol(1:2);
else           szvol = szvol(1:3);
end

nv = prod(szvol);

if(nv ~= size(mat,2) )
  fprintf('ERROR: fast_mat2vol: szvol inconsistent with mat2d\n');
  fprintf('  nv = %d, size(mat,2) = %d\n',nv,size(mat,2));
  return;
end
nframes = size(mat,1);

vol = reshape(transpose(mat), [szvol(:)' nframes]);

return;
