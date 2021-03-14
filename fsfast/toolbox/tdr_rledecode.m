function v = tdr_rledecode(rle)
% v = tdr_rledecode(rle)
%
% Decodes a run-length encoded sequence.
%
%


%
% tdr_rledecode.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

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
