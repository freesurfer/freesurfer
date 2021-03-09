function [f, mristruct] = fast_ldbslice(bstem,sliceno)
% [f mristruct] = fast_ldbslice(bstem, <sliceno>)
%
% sliceno is the zero-based slice number.
% If sliceno is not given or if it is < 0 or [], then
% the volume is loaded. 
%
% mristruct is mri info from the bhdr file. See fast_mri_struct
% and fast_ldbhdr.
%
% A single slice is rows-cols-frames
% The volume is     rows-cols-slices-frames
%
% See also fast_svbslice, fast_mri_struct, fast_ldbhdr.
%
%


%
% fast_ldbslice.m
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

f = [];
mristruct = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('USAGE: [f mristruct] = fast_ldbslice(bstem,sliceno)\n');
  return;
end

[nslices nrows ncols nframes endian bext hdrdat] = fmri_bvoldim(bstem);
if(nslices == 0)
  fprintf('ERROR with bvolume %s\n',bstem);
  return;
end

if(nargin == 1) sliceno = -1; end
if(isempty(sliceno)) sliceno = -1; end

if(sliceno >= nslices)
  fprintf('ERROR: requested slice %d exceeds number of slices %d\n',...
	  sliceno,nslices);
  return;
end

if(sliceno >= 0 & length(sliceno) == 1)
  % Read in a single slice %
  fname = sprintf('%s_%03d.%s',bstem,sliceno,bext);
  f = fmri_ldbfile(fname);
elseif(sliceno < 0)
  % Read in the volume %
  f = zeros(nrows,ncols,nslices,nframes);
  for s = 0:nslices-1
    tmp = fast_ldbslice(bstem,s);
    f(:,:,s+1,:) = squeeze(tmp);
  end
end

if(nargout == 2)
  bhdr = sprintf('%s.bhdr',bstem);
  mristruct = fast_ldbhdr(bhdr);
end

return;
