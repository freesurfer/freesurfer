function Ftrap = tdr_trapspect2d(wx,wy,width,center,ntrap)
%
% Ftrap = tdr_trapspect2d(wx,wy,<width>,<center>,<ntrap>)
%
% Sythesize the continuous spectrum of a 2D stair-step trapezoid at
% radial frequnecies wx, wy. 
%
% width  - [xwidth  ywidth],  default is 20 30
% center - [xcenter ycenter], default is 0 0
% ntrap  - number stair steps, default is 2
%
%


%
% tdr_trapspect2d.m
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

Ftrap = [];

if(nargin < 2 | nargin > 5)
  fprintf('Ftrap = tdr_trapspect2d(wx,wy,width,center,ntrap)\n');
  return;
end

if(exist('width') ~= 1) width = []; end
if(isempty(width)) width = [20 30]; end
if(length(width) ~= 2)
  fprintf('ERROR: width must have two components\n');
  return;
end

if(exist('center') ~= 1) center = []; end
if(isempty(center)) center = [0 0]; end
if(length(center) ~= 2)
  fprintf('ERROR: center must have two components\n');
  return;
end

if(exist('ntrap') ~= 1) ntrap = []; end
if(isempty(ntrap)) ntrap = 2; end

if(length(wx) ~= length(wy))
  fprintf('ERROR: wx and wy must be the same size\n');
  return;
end

% Synthsize the spectrum for a stair-step trapezoid in 1D
% For pure pulse, set xroll = 0;
x0 = center(1) - width(1)/2;;      % left edge
x1 = center(1) + width(1)/2;;      % right edge
if(ntrap > 0)
  xroll = [0:ntrap];
  nxroll = length(xroll);
  x00 = repmat(x0,[1 nxroll]);
  x00 = x00 + xroll;
  x11 = repmat(x1,[1 nxroll]);
  x11 = x11 - xroll;
else
  x00 = x0;
  x11 = x1;
end

Fxpulse = pulsespectrum(x00,x11,wx);
if(isempty(Fxpulse)) return; end
Fxtrap = mean(Fxpulse,2);

% Synthsize the spectrum for a stair-step trapezoid in 1D
% For pure pulse, set yroll = 0;
y0 = center(2) - width(2)/2;      % bottom edge
y1 = center(2) + width(2)/2;      % top edge
if(ntrap > 0)
  yroll = [0:ntrap];
  nyroll = length(yroll);
  y00 = repmat(y0,[1 nyroll]);
  y00 = y00 + yroll;
  y11 = repmat(y1,[1 nyroll]);
  y11 = y11 - yroll;
else
  y00 = y0;
  y11 = y1;
end

Fypulse = pulsespectrum(y00,y11,wy);
if(isempty(Fypulse)) return; end
Fytrap = mean(Fypulse,2);

Ftrap = Fxtrap .* Fytrap;

return;







