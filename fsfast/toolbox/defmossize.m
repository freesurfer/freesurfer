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
