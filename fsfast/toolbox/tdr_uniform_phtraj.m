function phi = tdr_uniform_phtraj(nksamples)
% phi = tdr_uniform_phtraj(nksamples)
%
% Assumes 'standard' siemens trajectory which:
%   Starts pos at +pi and goes neg. 
%   Ends at -pi+dphi
%   Phase decrements by dphi=2*pi/nksamples for each sample
%   Passes thru 0 at nksamples/2 + 1
%
%   If the direction is reversed, then you must just flipud(phi). 
%
% Note that this applies for both readout and phase encode.
%
%


%
% tdr_uniform_phtraj.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:35 $
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
  fprintf('phi = tdr_uniform_phtraj(nksamples)\n');
  return;
end

% Phase increment
dphi = 2*pi/nksamples;

% Starts at +pi, passes thru 0 at nc/2 + 1, ends at -pi+dphi
phi = dphi*[nksamples:-1:1]' - pi; 




