function cmap = nmrcolormap(ncmap,tail)
%
% cmap = nmrcolormap(ncmap,<tail>)
%


%
% nmrcolormap.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
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

if(nargin ~= 1 & nargin ~= 2)
  msg = 'Usage: cmap = nmrcolormap(ncmap,<tail>)';
  qoe(msg); error(msg);
end

if(nargin ~= 2) tail = 'pos'; end

ncmap2 = floor(ncmap/2);

z = [0:ncmap2-1]'/ncmap2; %'

r = zeros(ncmap,1);
g = zeros(ncmap,1);
b = zeros(ncmap,1);
g(ncmap2+1:ncmap) = z;

switch(tail)
  case {'pos','abs'},
    r(:) = 1;
    r(1:ncmap2) = z;
    cmap = [r g b];
  case {'neg'}
    b(:) = 1;
    b(1:ncmap2) = z;
    cmap = [r g b];
  case {'both','posneg','negpos'}
    cmappos = nmrcolormap(ncmap2,'pos');
    cmapneg = nmrcolormap(ncmap2,'neg');
    cmap = [cmapneg(ncmap2:-1:1,:); cmappos];
end


return;
