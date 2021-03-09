function [indTPExc, indTPInc] = fast_ldtpexcl(TPExcludeFile,TR,NTRs,nSkip)
% [indTPExc indTPInc] = fast_ldtpexcl(TPExcludeFile,TR,NTRs,nSkip)
%
% Returns a list of indices of time points which should
% be excluded and a list of those which should be included. These lists
% are based on what is in the Time Point Exlude file and on the value
% nSkip.  The TPExcludeFile variable can be empty.
%


%
% fast_ldtpexcl.m
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

if(nargin ~= 4)
  msg = 'USAGE: [indTPExc indTPInc] = fast_ldtpexcl(TPExcludeFile,TR,NTRs,nSkip)';
  qoe(msg);error(msg);
end

if(~isempty(TPExcludeFile))
  fid = fopen(TPExcludeFile);
  if(fid == -1)
    msg = sprintf('ERROR: could not open %s',TPExcludeFile);
    qoe(msg);error(msg);
  end
  TPX = fscanf(fid,'%f');
  fclose(fid);

  if(isempty(TPX)) 
    if(nSkip < 1)
      indTPExc = [];
      indTPInc = [1:NTRs];
    else
      indTPExc = [1:nSkip];
      indTPInc = [nSkip+1:NTRs];
    end
    return;
  end
  indTPExc = round(TPX/TR) + 1;
else
  indTPExc = [];
end

iOut = find(indTPExc < 1 | indTPExc > NTRs);
if(~isempty(iOut))
  msg = sprintf('TPX: %s has values out of bounds',TPExcludeFile);
  qoe(msg);error(msg);
end

if(nSkip > 0) indTPExc = unique([indTPExc; [1:nSkip]']); end %'

tmp = zeros(NTRs,1);
tmp(indTPExc) = 1;
indTPInc = find(tmp == 0);

return;






