function r = fmri_touch(fname);
% r = fmri_touch(fname);
% 
% Simple function to create a file called fname.  This is supposed
% to be something like the unix touch, but it has
% no effect if fname already exists.
%
%


%
% fmri_touch.m
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
  msg = 'USAGE: r = fmri_touch(fname);';
  qoe(msg); error(msg);
end

fname = deblank(fname);

fid = fopen(fname,'a');
if(fid == -1) 
  msg = sprintf('Could not open %s for appending',fname);
  qoe(msg); error(msg);
end

fclose(fid);

r = 0;
return;
