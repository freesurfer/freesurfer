function isMGH = MRIisMGH(fspec)
% isMGH = MRIisMGH(fspec)
%
% Determines whether the given file spec is an MGH/MGZ
% file based soley on its extension
%
% Returns non-zero if fspec is the name of an MGH file,
%   Returns 1 if it has a .mgh extension.
%   Returns 2 if it has a .mgz extension.
% Returns 0 otherwise.
%
% $Id: MRIisMGH.m,v 1.1 2004/11/09 19:14:57 greve Exp $

isMGH = [];
if(nargin ~= 1)
  fprintf('isMGH = MRIisMGH(fspec)\n');
  return;
end

isMGH = 0;
if(length(fspec) < 5) return; end

ext = fspec(end-3:end);

if(strcmp(ext,'.mgh')) isMGH = 1; end
if(strcmp(ext,'.mgz')) isMGH = 2; end

return;













