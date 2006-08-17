function [mcvhist, mcvunique] = bmm_mcvhist(mcvect)
% [mcvhist, mcvunique] = bmm_mcvhist(mcvect)
%
% Not tested yet!
% 
% Computes sparse histogram of the modified count vector (mcv) use
% in the two-method dependent form of the binomial mixture model.
%
% mcvect is the spatialdims-by-4 modified count vector as generated
% by bmm_mcvect.
%
% mcvunique - list of unique (sparse) vectors in mcvect
% mcvhist - number of times each unique vector appears in mcvect.
%   sum(mcvhist) = nvoxels
%
% Based on appendix of Genovese, et al, 1997. Estimating Test-Retest
% Reliability in Functional MR Imaging I: Statistical Methodology.
% MRM 38:497-507. j=mcvunique and nj=mcvhist.
%
% $Id: bmm_mcvhist.m,v 1.1 2006/08/17 05:21:52 greve Exp $
%

mcvhist=[];
mcv = [];
if(nargin ~= 1)
  fprintf('[mcvhist, mcvunique] = bmm_mcvhist(mcvect)\n');
  return;
end

sz = size(mcvect);
if(sz(end) ~= 4)
  fprintf('ERROR: mcvect must have 4 elements in last dim\n');
  return;
end

vsz = sz(1:end-1); % spatial dims
nv = prod(vsz);    % number of voxels

% Reshape to 2D nv-by-4
mcvect = reshape(mcvect,[nv 4]);

% Sort the rows in ascending order
mcvect = sortrows(mcvect);

% Get the unique vectors (nunique-by-4)
mcvunique = unique(mcvect,'rows');

% Now the tricky part. Get the number of times each of the unqiue
% vectors is replicated in mcvect. This is the histogram.
d = diff(mcvect);
a = sum(abs(d),2) > 0;
ind = [0; find(a); nv];
mcvhist = diff(ind);

return;
