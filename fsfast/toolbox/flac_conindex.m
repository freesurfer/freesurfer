function ind = flac_conindex(conname,flac)
% ind = flac_conindex(conname,flac)
% Returns the index of the given contrast name in the flac
% $Id: flac_conindex.m,v 1.1 2004/11/30 23:54:57 greve Exp $

ind = [];
if(nargin ~= 2)
  fprintf('ind = flac_conindex(conname,flac)\n');
  return;
end

ncon = length(flac.con);
for nthcon = 1:ncon
  if(strcmp(flac.con(nthcon).name,conname))
    ind = nthcon;
    return;
  end
end

% Will return empty ind if gets here.

return;
















