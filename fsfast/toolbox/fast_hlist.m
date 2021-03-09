function [avglist, stdlist] = fast_hlist(Nc,Nh)
% [avglist stdlist] = fast_hlist(Nc,Nh)


%
% fast_hlist.m
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

avglist = [];
stdlist = [];

if(nargin ~= 2)
  msg = '[avglist stdlist] = fast_hlist(Nc,Nh)';
  qoe(msg); error(msg);
end 

n = 1;
for c = 1:Nc
  for s = 1:2
    for h = 1:Nh
      if(s==1) avglist = [avglist n];
      else     stdlist = [stdlist n];
      end
      n = n + 1;
    end
  end
end

return;
