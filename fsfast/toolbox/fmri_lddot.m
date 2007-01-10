function dot = fmri_lddot(dotfile)
% dot = fmri_lddot(dotfile)


%
% fmri_lddot.m
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
  msg = 'USAGE: dot = fmri_lddot(dotfile)';
  qoe(msg);error(msg);
end

dot = load(dotfile);
[n1 n2] = size(dot);
ntp = n1/2;
nsensors = n2;

dot = reshape(dot, [ntp 2 nsensors]);
dot = permute(dot, [2 3 1]);

return;

