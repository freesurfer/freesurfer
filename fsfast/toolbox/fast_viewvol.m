function h = fast_viewvol(vol,volres,frameres)
% h = fast_viewvol(vol,volres,frameres)


%
% fast_viewvol.m
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

h = [];

if(nargin < 1 | nargin > 3)
  fprintf('USAGE: h = fast_viewvol(vol,volres,frameres)\n');
  return;
end

if(nargin < 2)      volres = [1 1 1]; end
if(isempty(volres)) volres = [1 1 1]; end

if(nargin < 3)        frameres = 1; end
if(isempty(frameres)) frameres = 1; end

