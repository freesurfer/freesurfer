function [avglist, stdlist] = fast_hlist(Nc,Nh)
% [avglist stdlist] = fast_hlist(Nc,Nh)


%
% fast_hlist.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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
