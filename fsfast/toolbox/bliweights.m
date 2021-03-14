function [ind, wind] = bliweights(szmtx, rc_hat)
% [ind wind] = bliweights(szmtx, rc_hat)
%
% Weights and indices for bilinear interpolation.
%
% szmtx is number of rows and columns in matrix
% rc_hat - off-grid rows and columns to interpolate (N x 2)
%
% ind (N x 4)
% w (N x 4)
%
% See also: uliweights, tliweights


%
% bliweights.m
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
  msg = 'USAGE: [ind wind] = bliweights(szimg, rc_hat)';
  error(msg);
end

szmtx = szmtx(1:2);
r_hat = rc_hat(:,1);
c_hat = rc_hat(:,2);
clear rc_hat;

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

% Compute the indices of each corner 
% of the interpolation square (row major!)
ind1 = sub2ind(szmtx,r1,c1);
ind2 = sub2ind(szmtx,r2,c1);
ind3 = sub2ind(szmtx,r1,c2);
ind4 = sub2ind(szmtx,r2,c2);
ind = [ind1 ind2 ind3 ind4];

% Compute the corresponding weights %
w1 = wr1 .* wc1;
w2 = wr2 .* wc1;
w3 = wr1 .* wc2;
w4 = wr2 .* wc2;

wind = [w1 w2 w3 w4];

return;



%------------------------------------------------------%
%------------------------------------------------------%
%------------------------------------------------------%
if(0) %---------------------------------------%
% Check for rows that are out of bounds %
roob = find(r_hat < 1 | r_hat  > szmtx(1));
if(~isempty(roob))
  msg = sprintf('Found %d row subscripts out of bounds',length(roob));
  error(msg);
end

% Check for columns that are out of bounds %
coob = find(c_hat < 1 | c_hat  > szmtx(2));
if(~isempty(coob))
  msg = sprintf('Found %d column subscripts out of bounds',length(coob));
  error(msg);
end

% Compute row upper and lower bounds %
r1 = floor(r_hat);
r2 = r1 + 1;
roob = find(r2 > szmtx(1));
r1(roob) = r1(roob) - 1;
r2(roob) = r2(roob) - 1;

% Compute column upper and lower bounds %
c1 = floor(c_hat);
c2 = c1 + 1;
coob = find(c2 > szmtx(2));
c1(coob) = c1(coob) - 1;
c2(coob) = c2(coob) - 1;

% Compute the corresponding weights %
u = (r_hat - r1);
t = (c_hat - c1);
w1 = (1-t).*(1-u);
w2 = (1-t).*u;
w3 = t .* (1-u);
w4 = t.*u;

end %---------------------------------------%
