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
%  The first  column is the sine   at the fundamental.
%  The second column is the cosine at the fundamental.
%
% Note: this really has nothing to do with the fast fourier
% transform.
%
% $Id: fast_fourier_reg.m,v 1.1 2004/06/03 17:53:51 greve Exp $
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
  X = [X sin(ph) cos(ph)];
end

return;





