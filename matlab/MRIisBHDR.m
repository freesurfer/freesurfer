function [isBHDR, fstem] = MRIisBHDR(fspec)
% [isBHDR fstem] = MRIisBHDR(fspec)
%
% Determines whether the given file spec is an BHDR
% file based soley on its extension. BHDR is a bshort or
% bfloat.
%
% Returns non-zero if fspec is the name of an BHDR file,
%   Returns 1 if it has a .bhdr extension.
% Returns 0 otherwise.
%
% $Id: MRIisBHDR.m,v 1.1 2004/11/13 16:47:36 greve Exp $

fstem = [];
isBHDR = [];
if(nargin ~= 1)
  fprintf('[isBHDR fstem] = MRIisBHDR(fspec)\n');
  return;
end

isBHDR = 0;
if(length(fspec) < 6) return; end

ext = fspec(end-4:end);

if(strcmp(ext,'.bhdr')) 
  isBHDR = 1; 
  fstem = fspec(1:end-5);
end

return;













