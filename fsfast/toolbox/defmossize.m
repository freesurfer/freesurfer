function tszmos = defmossize(N, tszmos)
% tszmos = defmossize(N, <tszmos>)
%
% Computes default mosaic size in terms of number of tile rows and
% colums.
%
% Inputs:
%   N - minimum number of slices from which to construct the mosaic.
%   tszmos - mosaic size [rows cols] in tiles (optional).
%
% Returns:
%   tszoms - mosaic size [rows cols] in tiles.
%
% The tszmos input argument can be used to specify the number of 
% tile rows or the number of tile colums.


%
% defmossize.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:29 $
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

if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: tszmos = defmossize(N, <tszmos>)';
  error(msg);
end

if(N == 0) 
  msg = 'Cannot make a mosaic with zero slices';
  error(msg);
end 

if(nargin == 1 | (nargin == 2 & isempty(tszmos)))
  Ntr = floor(sqrt(N));
  Ntc = ceil(N/Ntr);
  tszmos = [Ntr Ntc];
end

if(tszmos(1) == 0 & tszmos(2) ~= 0)  
  tszmos(1) = ceil(N/tszmos(2));
end

if(tszmos(1) ~= 0 & tszmos(2) == 0)  
  tszmos(2) = ceil(N/tszmos(1));
end

if(prod(tszmos) < N)
  msg = sprintf('Mosaic size (%d,%d) is not big enough for %d slices',...
                tszmos(1),tszmos(2),N);
  error(msg);
end


return;
