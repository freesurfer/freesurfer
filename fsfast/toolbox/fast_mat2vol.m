function vol = fast_mat2vol(mat,szvol,sliceflag)
%
% vol = fast_mat2vol(mat,szvol,<sliceflag>)
% 
% Reshapes a 2d matrix of size nframes X nvoxels
% to a volume of size [szvol nframes].
%
% szvol = [nrows ncols nslices]. If szvol only has two elements,
% then nslices=1. If sliceflag is set, then nslices=1, and the
% number of spatial voxels is assumed to be nrows*ncols.
%
% See also fast_vol2mat.
%
% $Id: fast_mat2vol.m,v 1.1 2004/04/28 18:41:39 greve Exp $

vol = [];

if(nargin < 2 | nargin > 3)
  fprintf('vol = fast_mat2vol(mat,szvol,<sliceflag>)\n');
  return;
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
