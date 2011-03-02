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
%


%
% bmm_mcvhist.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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
