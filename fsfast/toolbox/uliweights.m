function [ind, wind] = uliweights(szmtx, i_hat)
% [ind, wind] = uliweights(szmtx, i_hat)
%
% Weights and indices for (uni)linear interpolation from a 
% uniform grid.
%
% szmtx is number of elements in the grid.
% i_hat - off-grid indices to interpolate (N x 1)
%
% ind (N x 2)
% w (N x 2)
%
% See also: bliweights, tliweights


%
% uliweights.m
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

if(nargin ~= 2)
  msg = 'USAGE: [ind wind] = uliweights(szimg, i_hat)';
  error(msg);
end

szmtx = szmtx(1);

% Check for indices that are out of bounds %
ioob = find(i_hat < 1 | i_hat  > szmtx(1));
if(~isempty(ioob))
  msg = sprintf('Found %d subscripts out of bounds',length(ioob));
  error(msg);
end

% Compute upper and lower on-grid indices %
i1 = floor(i_hat);
i2 = i1 + 1;
ioob = find(i2 > szmtx(1));
i1(ioob) = i1(ioob) - 1;
i2(ioob) = i2(ioob) - 1;
ind = [i1 i2];

% Compute the corresponding weights %
w2 = i_hat - i1;
w1 = 1 - w2;
wind = [w1 w2];

return;
