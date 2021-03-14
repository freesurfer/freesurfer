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
%


%
% fast_svbslice.m
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

err = 1;

if(nargin < 2 | nargin > 5)
  fprintf('err = fast_svbslice(y,stem,sliceno,<bext>,<mristruct>)\n');
  return;
end

if(exist('sliceno') ~= 1) sliceno = -1; end
if(exist('bext') ~= 1) bext = ''; end
if(isempty(bext)) bext = 'bfloat'; end

if(strcmp(bext,'bfloat') == 0 & strcmp(bext,'bshort') == 0)
  fprintf('ERROR: bext = %s, must be bfloat or bshort\n',bext);
  return;
end
  
if(sliceno >= 0)
  % Save as a single slice %
  fname = sprintf('%s_%03d.%s',stem,sliceno,bext);
  err = fmri_svbfile(y,fname);
  if(err) return; end
else
  % Save as an entire volume %
  nrows = size(y,1);
  ncols = size(y,2);
  nslices = size(y,3);
  nframes = size(y,4);
  yslice = zeros(nrows,ncols,nframes);
  for slice = 0:nslices-1
    fname = sprintf('%s_%03d.%s',stem,slice,bext);
    % Cant squeeze y in case nrows or ncols == 1
    yslice = y(:,:,slice+1,:);
    err = fmri_svbfile(yslice,fname);
    if(err) return; end
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
if(exist('mristruct') == 1 & ~isempty(mristruct))
  szf = size(y);
  fdims = length(szf);
  if(sliceno < 0) % volume
    if(fdims == 3) nframes = 1;
    else           nframes = szf(4);
    end
    mristruct.voldim = szf([2 1 3]); % cols, rows, slices
  else % slice
    if(fdims == 2) nframes = 1;
    else           nframes = szf(3);
    end
    mristruct.voldim(1:2) = szf([2 1]); % cols, rows
  end
  mristruct.nframes = nframes;
  bhdrfile = sprintf('%s.bhdr',stem);
  fast_svbhdr(mristruct,bhdrfile);
end

err = 0;

return;







