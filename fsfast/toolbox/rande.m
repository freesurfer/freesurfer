function r = rande(rdim)
%
% r = rande(dim)
%
% Generates matrix of given dimension whose elements 
% are exponentially distributed. The mean and std over 
% all the elements asymptotically approach one.
%
% $Id: rande.m,v 1.1 2003/03/04 20:47:41 greve Exp $

if(nargin == 0) rdim = 1; end

r = -log(rand(prod(rdim),1));
r = reshape(r, [rdim 1]);

return;
