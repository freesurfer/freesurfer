function v = tdr_rledecode(rle)
% v = tdr_rledecode(rle)
%
% Decodes a run-length encoded sequence.
%
% $Id: tdr_rledecode.m,v 1.1 2006/05/26 23:49:11 greve Exp $

nrle = length(rle);

v = rle(1);
nthrle = 2;
while(nthrle <= nrle)

  val = rle(nthrle);
  if(rle(nthrle) == rle(nthrle-1))
    % repeated value
    nthrle = nthrle + 1;
    nreps = rle(nthrle);
    nthrle = nthrle + 1;
    v = [v repmat(val,[1 nreps+1])];
  else
    % not a repeated value
    nthrle = nthrle + 1;
    v = [v val];
  end
  
end

return;
