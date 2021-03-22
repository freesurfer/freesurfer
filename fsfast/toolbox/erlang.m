function y = erlang(x,avg,r)
% y = erlang(x,avg,r)
%
% Generates theoretical erlang distribution with average avg and order
% r. If avg and r are not specfied, both are assumed to be 1. The
% variance will be (avg^2)/r. The stddev = avg/sqrt(r).
%
% To test emperically:
%   y = rande(100000,5) + 2.3; % avg is 2.3+1 = 3.3
%   [h x] = hist(y,100);  h = h/max(h);
%   hest = erlang(x,2.3+1,5); hest = hest/max(hest);
%   plot(x,h,x,hest);
%
% y = r*((r*(x-(avg-1)))^(r-1)) * exp(-r*(x-(avg-1))) / (r-1)!
%
% See also rande.
%
%


%
% erlang.m
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

y = [];

if(nargin < 1 | nargin > 3)
  fprintf('y = erlang(x,avg,r)\n');
  return;
end

if(exist('avg')~=1) avg = 1; end
if(exist('r')~=1)  r = 1; end

x = x - (avg-1);
y = zeros(size(x));
indgez = find(x >= 0);

y(indgez) = r*((r*x(indgez)).^(r-1)) .* exp(-r*x(indgez)) / factorial(r-1);

return;
  
  
  
  
