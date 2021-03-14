function [havg, hstd, Nh] = hsa_convert(hsa,Nnnc)


%
% hsa_convert.m
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

