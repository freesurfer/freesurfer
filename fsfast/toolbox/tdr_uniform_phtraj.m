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

if(nargin ~= 1)
  fprintf('phi = tdr_uniform_phtraj(nksamples)\n');
  return;
end

% Phase increment
dphi = 2*pi/nksamples;

% Starts at +pi, passes thru 0 at nc/2 + 1, ends at -pi+dphi
phi = dphi*[nksamples:-1:1]' - pi; 




