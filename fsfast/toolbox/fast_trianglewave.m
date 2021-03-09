function [w period] = fast_trianglewave(t,spec,spectype)
% [w period] = fast_trianglewave(t,period)
% [w period] = fast_trianglewave(t,period,'period')
% [w period] = fast_trianglewave(t,slope,'slope')
%
% Produces a triangle wave of the given period. Starts at 0 at t=0
% peaks at 1 at t=period/2, then returns to 0 at t=period.
% If the slope is specified, then period = 2/abs(slope). If 
% slope=0, then period=inf and w will be set to 0.

%
% fast_trianglewave.m
%
% Original Author: Douglas Greve
% Jan 29, 2019
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

if(nargin < 2 | nargin > 3)
  fprintf('w = fast_trianglewave(t,period)\n');
  fprintf('w = fast_trianglewave(t,period,''period'')\n');
  fprintf('[w period] = fast_trianglewave(t,slope,''slope'')\n')
  return;
end

if(nargin == 2)
  period = spec;
  slope = 1;
else
  if(strcmp(spectype,'period'))
    period = spec;
    slope = 1;
  else
    slope = spec;
    if(slope==0)
      period = inf;
      w = zeros(size(t));
      return;
    end
    % 1/slope is time to get to a value of 1
    period = 2/abs(slope);
  end
end

% This is a number between 0 and 1
w = rem(t,period)/period;
% Take all values greater than 0.5 and make it 1-value
% so that all values are between 0 and 0.5
ind = find(w>0.5);
w(ind) = 1-w(ind);
% Rescale so that all values are between 0 and 1
w = (2*sign(slope))*w;

return








