function err = fast_svbslice(y,stem,sliceno,bext,mristruct)
% err = fast_svbslice(y,stem,sliceno,bext,mristruct)
%
% size(y) = [rows cols frames]
% sliceno is zero-based
% if bext is not specificed or is null, it is set to bfloat
% mristruct (see fast_mri_struct) is used to create stem.bhdr
%
% if sliceno is < 0 then size(y) = [rows cols slices frames]
% and each slice is saved.
%
% See also fast_ldbslice, fast_mri_struct, fast_svbhdr.
% 
% $Id: fast_svbslice.m,v 1.6 2003/08/04 19:13:31 greve Exp $

err = 1;

if(nargin < 2 | nargin > 5)
  fprintf('err = fast_svbslice(y,stem,sliceno,<bext>,<mristruct>)\n');
  return;
end

if(exist('bext') ~= 1) bext = ''; end
if(isempty(bext)) bext = 'bfloat'; end

if(strcmp(bext,'bfloat') == 0 & strcmp(bext,'bshort') == 0)
  fprintf('ERROR: bext = %s, must be bfloat or bshort\n',bext);
  return;
end
  
if(sliceno >= 0)
  % Save as a single slice %
  fname = sprintf('%s_%03d.%s',stem,sliceno,bext);
  fmri_svbfile(y,fname);
else
  % Save as an entire volume %
  nslices = size(y,3);
  for slice = 0:nslices-1
    fname = sprintf('%s_%03d.%s',stem,slice,bext);
    fmri_svbfile(squeeze(y(:,:,slice+1,:)),fname);
  end

  % Check for extra slices and delete them %
  slice = nslices;
  while(1)
    fname = sprintf('%s_%03d.%s',stem,slice,bext);
    fid = fopen(fname,'r');
    if(fid ~= -1) delete(fname);
    else break;
    end
    fname = sprintf('%s_%03d.hdr',stem,slice);
    fid = fopen(fname,'r');
    if(fid ~= -1) delete(fname);
    else break;
    end
    slice = slice+1;
  end

end

% Save the bhdr file %
if(exist('mristruct') == 1)
  szf = size(y);
  fdims = length(szf);
  nframes = szf(fdims);
  mristruct.nframes = nframes;
  bhdrfile = sprintf('%s.bhdr',stem);
  fast_svbhdr(mristruct,bhdrfile);
end

err = 0;

return;







