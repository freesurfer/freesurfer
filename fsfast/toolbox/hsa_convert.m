function [havg, hstd, Nh] = hsa_convert(hsa,Nnnc)


%
% hsa_convert.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.3 $
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



  Nc = Nnnc+1;

  if(size(hsa,2) == 1 & length(size(hsa)) == 2)
    ud.nRows = 1;
    ud.nCols = 1;
    Nch2 = size(hsa,1);
  else
    [ud.nRows ud.nCols Nch2] = size(hsa);
  end
  Nh = Nch2/(2*Nc);

  hsa2 = permute(hsa, [3 2 1 4]); % t,c,r
  hsa3 = reshape(hsa2, [Nh 2 Nc ud.nCols ud.nRows ]); % h stat cond col row
  clear hsa2;

  hsa4 = permute(hsa3, [5 4 3  1 2 ]); % row col cond h stat
  clear hsa3;

  havg = (hsa4(:,:,:,:,1)); % col row cond havg 
  hstd = (hsa4(:,:,:,:,2)); % col row cond hstd

  %havg = squeeze(hsa4(:,:,:,:,1)); % col row cond havg 
  %hstd = squeeze(hsa4(:,:,:,:,2)); % col row cond hstd
  clear hsa4;

return

