function [slc, slcdel] = fmri_ldslicedelay(fname)
% [slc slcdel] = 


%
% fmri_ldslicedelay.m
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

if(nargin ~= 1)
  msg = 'USAGE: [slc slcdel] = fmri_ldslicedelay(fname)';
  qoe(msg);error(msg);
end

fname = deblank(fname);

fid = fopen(fname,'r');
if(fid == -1)
  msg = sprintf('Could not open %s',fname);
  qoe(msg);error(msg);
end

tmp = fscanf(fid,'%f');
fclose(fid);

nslices = length(tmp)/2;
tmp = reshape(tmp, [2 nslices])'; %'
slc    = tmp(:,1);
slcdel = tmp(:,2);

%% resort by slice number %%
[y i] = sort(slc);
slc    = y;
slcdel = slcdel(i);

return;
