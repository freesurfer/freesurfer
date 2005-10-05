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
% $Id: fast_mat2vol.m,v 1.4 2005/10/05 20:43:55 greve Exp $

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
    fprintf('ERROR: need szvol if mat is not an mri struct\n');
    return;
  end
end

if(length(szvol) < 2)
  fprintf('ERROR: szvol dim must be >= 2\n');
  return;
end

if(~exist('sliceflag','var')) sliceflag = []; end
if(isempty(sliceflag)) sliceflag = 0; end

if(sliceflag)  szvol = szvol(1:2);
else           szvol = szvol(1:3);
end

nv = prod(szvol);

if(nv ~= size(mat,2) )
  fprintf('ERROR: szvol inconsistent with mat2d\n');
  return;
end
nframes = size(mat,1);

vol = reshape(mat', [szvol nframes]);

return;
