function kimg2 = tdr_kshift(kimg,shift,dim)
% kimg2 = tdr_kshift(kimg,shift,<dim>)
%
% Shift in space of shift voxels (can be fractional),
% implemented by applying a phase roll to the kspace
% image. kimg can have multiple slices and/or frames.
%
% dim is the dimension in which to apply the shift
%   dim = 1 -- apply along the columns (default)
%   dim = 2 -- apply along the rows
%
% To shift in the phase encoding direction, dim = 1
%
% I'm not sure whether dim=2 actually works. 10/27/03
%
%
%


%
% tdr_kshift.m
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

kimg2 = [];
if(nargin ~= 2 & nargin ~= 3)
  fprintf('kimg2 = tdr_kshift(kimg,shift,dim)\n');
  return;
end

if(exist('dim') ~= 1) dim = 1; end
if(dim ~= 1 & dim ~= 2)
  fprintf('ERROR: dim = %d, must be 1 or 2\n',dim);
  return;
end

if(dim > length(size(kimg)))
  fprintf('ERROR: dim = %d, exceeds matrix dim\n',dim);
  return;
end

[nrows ncols nslices nframes] = size(kimg);
szkimg = [nrows ncols nslices nframes];

N = szkimg(dim);
ph = 2*pi * [0:N-1]'/N;
ph = ph - ph(N/2+1); 
vdelta = exp(-i*ph*shift);

if(dim == 1)
  kimg2 = kimg .* repmat(vdelta, [1  ncols nslices nframes]);
else
  kimg2 = kimg .* repmat(vdelta', [nrows 1 nslices nframes]);
end

return;








