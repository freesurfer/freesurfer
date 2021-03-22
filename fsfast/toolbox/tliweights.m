function [ind, wind] = tliweights(szmtx, rcs_hat)
% [ind wind] = tliweights(szmtx, rcs_hat)
%
% Weights and indices for trilinear interpolation.
%
% szmtx is number of rows, columns, and slices in volume (one-based)
% rcs_hat - off-grid row, column, and slice to interpolate (N x 3)
%
% ind (N x 8) -- on-grid indices (row major, one-based)
% w (N x 8) -- corresponding weights 
%
% See also: uliweights, bliweights


%
% tliweights.m
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
  msg = 'USAGE: [ind wind] = tliweights(szimg, rcs_hat)';
  error(msg);
end

szmtx = szmtx(1:3);
r_hat = rcs_hat(:,1);
c_hat = rcs_hat(:,2);
s_hat = rcs_hat(:,3);
clear rcs_hat;

% Get the upper and lower on-grid row subscripts %
% and their weights from unilinear interpolation. %
[r_ind wr_ind] = uliweights(szmtx(1), r_hat);
r1 = r_ind(:,1);
r2 = r_ind(:,2);
wr1 = wr_ind(:,1);
wr2 = wr_ind(:,2);
clear r_ind wr_ind;

% Get the upper and lower on-grid column subscripts %
% and their weights from unilinear interpolation. %
[c_ind wc_ind] = uliweights(szmtx(2), c_hat);
c1 = c_ind(:,1);
c2 = c_ind(:,2);
wc1 = wc_ind(:,1);
wc2 = wc_ind(:,2);
clear c_ind wc_ind;

% Get the upper and lower on-grid slice subscripts %
% and their weights from unilinear interpolation. %
[s_ind ws_ind] = uliweights(szmtx(3), s_hat);
s1 = s_ind(:,1);
s2 = s_ind(:,2);
ws1 = ws_ind(:,1);
ws2 = ws_ind(:,2);
clear s_ind ws_ind;

% Compute the indices of each corner 
% of the interpolation cube (row major!)
ind1 = sub2ind(szmtx,r1,c1,s1);
ind2 = sub2ind(szmtx,r2,c1,s1);
ind3 = sub2ind(szmtx,r1,c2,s1);
ind4 = sub2ind(szmtx,r2,c2,s1);
ind5 = sub2ind(szmtx,r1,c1,s2);
ind6 = sub2ind(szmtx,r2,c1,s2);
ind7 = sub2ind(szmtx,r1,c2,s2);
ind8 = sub2ind(szmtx,r2,c2,s2);
ind = [ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8];

% Compute the corresponding weights %
w1 = wr1 .* wc1 .* ws1;
w2 = wr2 .* wc1 .* ws1;
w3 = wr1 .* wc2 .* ws1;
w4 = wr2 .* wc2 .* ws1;
w5 = wr1 .* wc1 .* ws2;
w6 = wr2 .* wc1 .* ws2;
w7 = wr1 .* wc2 .* ws2;
w8 = wr2 .* wc2 .* ws2;

wind = [w1 w2 w3 w4 w5 w6 w7 w8];

return;
