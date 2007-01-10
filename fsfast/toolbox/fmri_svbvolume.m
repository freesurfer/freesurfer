function fmri_ldbvolume(vol,stem,voldim,ext)
% fmri_ldbvolume(vol,stem,voldim,ext)
%
% Saves a volume in bfile format where vol is a 4D structure
% of dimension Nslices X Nrows X Ncols X Ndepth.
%
% '$Id: fmri_svbvolume.m,v 1.2 2007/01/10 22:02:33 nicks Exp $'


%
% fmri_svbvolume.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
%    $Revision: 1.2 $
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

if(nargin ~= 2 & nargin ~= 3 & nargin ~= 4)    
  fprintf(2,'USAGE: fmri_svbvolume(vol,stem,<voldim>,<ext>)\n');
  qoe; error; return;
end

if(nargin > 2)  
  if(~isempty(voldim))  
    vol = reshape(vol, voldim); 
  else
    voldim = size(vol);
  end
else
  voldim = size(vol);
end

if(nargin ~= 4)  ext = 'bfloat'; end

nslices = voldim(1);
for slice = 0:nslices-1
  % fprintf('Saving Slice %3d\n',slice);
  fname = sprintf('%s_%03d.%s',stem,slice,ext);
  z = shiftdim(vol(slice+1,:,:,:),1);
  fmri_svbfile(z, fname);
end

return;
