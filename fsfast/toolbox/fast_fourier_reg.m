function X = fast_fourier_reg(period,nf,TR,nharmonics)
% X = fast_fourier_reg(period,nf,<TR>,<nharmonics>)
%
% Creates a fourier regression matrix.
%
% period - stimulation period in seconds (see TR)
% nf - number of frames
% TR - TR in seconds. If not present or empty, TR=1s
% nharmonics - number of harmonics to add. If not present
%   or empty, nharmonics=0.
%
% X - matrix with nf rows and 2*(1+nharmonics) columns.
%  The first  column is the cosine (real) at the fundamental.
%  The second column is the sine   (imag) at the fundamental.
%
% Note: this really has nothing to do with the fast fourier
% transform.
%


%
% fast_fourier_reg.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2010/04/17 18:52:22 $
%    $Revision: 1.3 $
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

if(nargin < 2 | nargin > 4)
  fprintf('X = fast_fourier_reg(period,nf,<TR>,<nharmonics>)\n');
  return;
end

if(~exist('TR','var')) TR = []; end
if(isempty(TR)) TR = 1; end

if(~exist('nharmonics','var')) nharmonics = []; end
if(isempty(nharmonics)) nharmonics = 0; end

t = TR*[0:nf-1]';

X = [];
for n = 0:nharmonics
  f = (n+1)/period;
  ph = 2*pi*t*f;
  X = [X cos(ph) sin(ph)]; % Real then Imag
end

return;





