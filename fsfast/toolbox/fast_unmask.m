function vol = fast_unmask(data,imask,volsize)
% vol = fast_unmask(data,imask,volsize)
%
% data is nf-by-nmask
% imask has nmask elements
% Creates vol which is ns-nr-nc-nf with zeros
%  in places without mask. 
% If imask is empty, converts data into a volume.
%


%
% fast_unmask.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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

vol = [];

if(nargin ~= 3)
  fprintf('vol = fast_unmask(data,imask,volsize)\n');
  return;
end

if(length(volsize) < 3)
  fprintf('ERROR: volsize must have at least 3 dimensions \n');
  return;
end
volsize = volsize(1:3);
nv = prod(volsize);
[nf ndata] = size(data);

if(isempty(imask))
  if(ndata ~= nv)
    fprintf('ERROR: ndata ~= nv\n');
    return;
  end
  vol = reshape(data', [volsize nf]); %'
  return;
end

if(max(imask) > nv)
  fprintf('ERROR: max(imask) > nv \n');
  return;
end

if(size(data,2)==1) data = data'; end %'

nmask = length(imask);
if(ndata ~= nmask)
  fprintf('ERROR: ndata ~= nmask \n');
  return;
end

vol = zeros([nf nv]);
vol(:,imask) = data;
vol = reshape(vol', [volsize nf]); %'

return;
