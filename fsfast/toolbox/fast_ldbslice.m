function [f, nv, nf] = fast_ldbslice(bstem,sliceno)
% [f nv nf] = fast_ldbslice(bstem, <sliceno>)
% If sliceno is not given or if it is < 0, then
% the volume is loaded. nv is the number of spatial
% voxels (ie, nrows*cols if a slice or nrows*cols*nslice 
% if a volume). nf is the number of frames.

f = [];
nv = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('USAGE: [f nv nf] = fast_ldbslice(bstem,sliceno)\n');
  return;
end

[nslices nrows ncols nt endian bext hdrdat] = fmri_bvoldim(bstem);
if(nslices == 0)
  fprintf('ERROR with bvolume %s\n',bstem);
  return;
end

if(nargin == 1) sliceno = -1; end

if(length(find(sliceno >= nslices))~=0)
  fprintf('ERROR: requested slice %d exceeds number of slices %d\n',...
	  sliceno(end),nslices);
  return;
end

if(sliceno(1) >= 0 & length(sliceno) == 1)
  fname = sprintf('%s_%03d.%s',bstem,sliceno,bext);
  f = fmri_ldbfile(fname);
  sz = size(f);
  nv = sz(1)*sz(2);
  if(length(sz)>2) nf = sz(3);
  else nf = 1;
  end
elseif(length(sliceno) > 1)
  nth = 1;
  for s = sliceno
    tmp = fast_ldbslice(bstem,s);
    f(nth,:,:,:) = squeeze(tmp);
    nth = nth+1;
  end
else
  f = fmri_ldbvolume(bstem);
  sz = size(f);
  nv = sz(1)*sz(2)*sz(3);
  if(length(sz)>3) nf = sz(4);
  else nf = 1;
  end
end

return;
