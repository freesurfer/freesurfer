function [x, a, b, c, ya, yb, yc] = gr1dmin_update(a, b, c, ya, yb, yc, x, yx)
% [x, a, b, c, ya, yb, yc] = gr1dmin_update(a, b, c, ya, yb, yc, x, yx)
%
% Computes updated parameters for Golden Ratio 1D minimization.
%
% Computes a new triplet and test point given the old triplet,
% test point, and value of the function at the test point based
% on the Golden Search 1D minimization (based on Numerical Recipces
% in C, p 398).
%
% a < b < c, ya < yb < yc
% a < x < c, ya < yx < yc
%
% This function can be invoked in one of two ways:
%
% (1) To get the next test point from a triplet (good for init)
%        x = gr1dmin_update(a, b, c) 
% (2) Get next test point and update triplets
%       [x a b c ya yb yc] = gr1dmin_update(a, b, c, ya, yb, yc, x, yx)
%
%


%
% gr1dmin_update.m
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

if(nargin ~= 3 & nargin ~= 8)
  fprintf(['ERROR: [x, a, b, c, ya, yb, yc] = gr1dmin_update(a, b, c, ya,' ...
	   ' yb, yc, x, yx)\n']);
  return;
end

if(nargin == 3 & nargout ~= 1)
  fprintf('ERROR: only one output with 3 inputs\n');
  return;
end

GoldenRatio = .61803399;
CompGoldenRatio = 1 - GoldenRatio;

if(nargin > 3)
  if( a < x & x < b)
    % x is between a and b %
    if(yx < yb)
      c  =  b;
      yc = yb;
      b  =  x;
      yb = yx;
    else
      a  =  x;
      ya = yx;
    end
  else
    % x is between b and c %
    if(yb < yx)
      c  =  x;
      yc = yx;
    else
      a  =  b;
      ya = yb;
      b  =  x;
      yb = yx;
    end
  end
end

if(abs(c-b) > abs(b-a))
  x = b + CompGoldenRatio*(c-b);
else
  x = b - CompGoldenRatio*(b-a);
end

return;



