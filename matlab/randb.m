function b = randb(p,ntrials,szb)
% b = randb(p,ntrials,<szb>)
% 
% Draw samples from a binomial distribution.
%
% p is the probability of a trial being "successful" 
%   (eg, probability that a voxel will be declared 
%   active (whether it is or not)). Scalar.
% ntrials - number of trials. Scalar.
% szb - size of b. Can be a vector. Default is 1.
%
% b will be a number between 0 and ntrials indicating the 
%   number of "successful" trials.
%
% Example:
%   b = randb(p,ntrials,100000);
%   x = [0:ntrials];
%   h = hist(b,x);
%   fx = binomialpdf(x,ntrials,p);
%   plot(x,fx,x,h/sum(h));
%


%
% randb.m
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


b = [];
if(nargin < 2 | nargin > 3)
  fprintf('b = randb(p,ntrials,szb)\n');
  return;
end

if(~exist('szb','var')) szb = []; end
if(isempty(szb)) szb = 1; end

szb = szb(:)';
nv = prod(szb);
b = reshape(sum(rand(ntrials,nv)<p),[szb 1]);

return;








