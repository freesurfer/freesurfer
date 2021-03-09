function C = mkprestimcon(sxa,wcond)
% C = mkprestimcon(sxa,wcond)


%
% mkprestimcon.m
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

C = [];
if(nargin ~= 2)
  fprintf('C = mkprestimcon(sxa,wcond)\n');
  return;
end

nPre = round(sxa.TPreStim/sxa.TER);
nFIR = sxa.Nh;
Nc = sxa.Nc-1;
wcond = wcond(:)';

if(length(wcond) ~= Nc)
  fprintf('ERROR: wcond has wrong number of items\n');
  return;
end

a = zeros(nFIR);
a(:,1:nPre) = -1/nPre;
b = a + eye(nFIR);
C = kron(wcond,b); 

return;
