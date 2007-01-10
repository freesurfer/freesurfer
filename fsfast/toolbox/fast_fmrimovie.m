function fmrimovie = fast_fmrimovie(f)
% fmrimovie = fast_fmrimovie(f)
%


%
% fast_fmrimovie.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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

if(nargin ~= 1)
  msg = 'USAGE: fmrimovie = fast_fmrimovie(f)';
  qoe(msg);error(msg);
end

figure;
imagesc(f(:,:,1));
axis image;
set(gca,'nextplot','replacechildren');

%f = diff(f,[],3);

% Adjust the scale to 1 - 64 %
fmax = max(reshape1d(f));
fmin = min(reshape1d(f));

f = 64*(f-fmin)/(fmax-fmin) + 1;

colormap(gray);

for j = 1:size(f,3);
  image(f(:,:,j));
  %colorbar;
  h = text(10,10,num2str(j));
  set(h,'color',[1 1 1])
  set(h,'fontsize',15)
  fmrimovie(j)=getframe;
end  

close

return
