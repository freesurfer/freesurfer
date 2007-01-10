function [Mss, Rss] = fmri_subsampmat(TR,Ntp,TS)
% [Mss Rss] = fmri_subsampmat(TR,Ntp,TS)
%
% Creates a subsampling matrix


%
% fmri_subsampmat.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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


if(nargin ~= 3)
  msg = 'USAGE: [Mss Rss] = fmri_subsampmat(TR,Ntp,TS)';
  qoe(msg);error(msg);
end

Rss = floor(TR/TS);
if(Rss < 1)
  msg = sprintf('Subsampling factor must be >= 1');
  qoe(msg);error(msg);
end

if(Rss == 1)
  Mss = eye(Ntp);
  return;
end

Nsp = Rss*Ntp;

Mss = zeros(Ntp,Nsp);
w   = [1:Rss:Nsp];
z   = [1:Ntp];
ind = sub2ind([Ntp Nsp], z,w);
Mss(ind) = 1;


return;
