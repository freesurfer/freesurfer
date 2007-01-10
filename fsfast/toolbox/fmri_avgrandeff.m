function [havg, hstd] = fmri_avgrandeff(hsum, hsum2, dof)
%
% [havg hstd] = fmri_avgrandeff(hsum, hsum2, dof)
%
% Comptutes the average and std dev for rand effects accumulated
% values.
%
%


%
% fmri_avgrandeff.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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

if(nargin ~= 3)
  msg = 'USAGE: [havg hstd] = fmri_avgrandeff(hsum, hsum2, dof)'
  qoe(msg);error(msg);
end

if(dof < 2)
  msg = sprintf('dof = %d, must be > 1',dof);
  qoe(msg);error(msg);
end

havg = hsum/dof;
hstd = real(sqrt( (hsum2 - dof*(havg.^2))/(dof-1)));

%tmp = (hsum2 - dof*(havg.^2))/(dof-1);
%keyboard
return;
