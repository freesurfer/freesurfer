function p = bvbinomialpdf(c,ntrials,p00,p01,p10)
% p = bvbinomialpdf(c,ntrials,p00,p01,p10)
%
% Not tested yet!
% 
% PDF of bivariate binomial distribution. Each trial consists of
% two coin flips, so there are four possible outcomes: 00, 01, 10,
% 11. The probabilities of only the three are needed because
% they must sum to 1. c is a 3-tuple row vector of the number of
% times the first three outcomes were realized over ntrials. Again
% the number fourth outcomes is not needed because the number of
% outcomes must sum to ntrials. 
%
% c - row vector of number of outcomes. Only the first 3 columns
% are used. If more than one row, then a p is computed for each
% row. 
%
% p - probability of seeing the distribution of c outcomes.
%
%


%
% bvbinomialpdf.m
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

p = [];
if(nargin ~= 5)
  fprintf('p = bvbinomialpdf(c,ntrials,p00,p01,p10)\n');
  return;
end

if(size(c,2) < 3)
  fprintf('ERROR: c must have at least 3 cols\n');
  return;
end

% Compute the prop of the 4th outcome
p11 = 1 - (p00 + p01 + p10);

% Turn c into a 4-tuple.
c = c(:,1:3);
csum = sum(c,2);
c = [c ntrials-csum];

nrows = size(c,1);
pvect = [p00 p01 p10 p11];
pp = repmat(pvect,[nrows 1]).^c;
p = prod(pp,2);

keyboard

return;

