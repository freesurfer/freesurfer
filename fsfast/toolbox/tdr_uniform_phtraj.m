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
% $Id: tdr_uniform_phtraj.m,v 1.1 2006/05/26 23:49:11 greve Exp $

if(nargin ~= 1)
  fprintf('phi = tdr_uniform_phtraj(nksamples)\n');
  return;
end

% Phase increment
dphi = 2*pi/nksamples;

% Starts at +pi, passes thru 0 at nc/2 + 1, ends at -pi+dphi
phi = dphi*[nksamples:-1:1]' - pi; 




