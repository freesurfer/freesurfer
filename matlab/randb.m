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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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








