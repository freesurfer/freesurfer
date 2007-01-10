function rimgstack = tdr_uniform_idft2d(kimgstack,revcol,ph0col,revrow,ph0row)
% rimgstack = tdr_uniform_idft2d(kimgstack,<revcol>,<ph0col>,<revrow>,<ph0row>)


%
% tdr_uniform_idft2d.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:35 $
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

if(nargin < 1 | nargin > 5)
  fprintf('rimgstack = tdr_uniform_idft2d(kimgstack,<revcol>,<ph0col>,<revrow>,<ph0row>)');
  return;
end

if(~exist('revcol','var')) revcol = []; end
if(isempty(revcol))        revcol = 0; end

if(~exist('ph0col','var')) ph0col = []; end
if(isempty(ph0col))        ph0col = 0; end

if(~exist('revrow','var')) revrow = []; end
if(isempty(revrow))        revrow = 0; end

if(~exist('ph0row','var')) ph0row = []; end
if(isempty(ph0row))        ph0row = 0; end

sz = size(kimgstack);
ndim = length(sz);
nrows = sz(1);
ncols = sz(2);
if(ndim > 2)
  nother = prod(sz(3:end));
else
  nother = 1;
end

% Encodes a column vector of rimg
Fcol = tdr_uniform_dft1d(nrows,nrows,revrow,ph0row); 
% Recons a column vector of kimg
Rcol = Fcol'/nrows; % = inv(Fcol'*Fcol)*Fcol'

% Encodes a row vector of rimg
Frow = tdr_uniform_dft1d(ncols,ncols,revcol,ph0col); 
% Recons a column vector of kimg
Rrow = Frow'/ncols; % = inv(Frow'*Frow)*Frow'

% Recon the column vectors first
kimgstack = reshape(kimgstack,[nrows ncols*nother]);
rimgstack = Rcol*kimgstack;
rimgstack = reshape(rimgstack,[nrows ncols nother]);

% Now recon the row vectors
rimgstack = permute(rimgstack,[2 1 3]);
rimgstack = reshape(rimgstack,[ncols nrows*nother]);
rimgstack = Rrow*rimgstack;
rimgstack = reshape(rimgstack,[ncols nrows nother]);
rimgstack = permute(rimgstack,[2 1 3]);

rimgstack = reshape(rimgstack,sz);

return;

