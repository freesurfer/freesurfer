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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:35 $
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
