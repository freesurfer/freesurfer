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
% $Id: fast_vol2mat.m,v 1.2 2005/09/15 15:53:45 greve Exp $

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

mat = reshape(vol,[nv nf])';


return;
