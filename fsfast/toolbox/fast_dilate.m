function md = fast_dilate(m,ndilations,erodeflag,flag2d)
% md = fast_dilate(m,<ndilations>,<erodeflag>,<flag2d>)
%
% Dilate or erode a binary mask.
%  m - input binary mask
%  ndilations - number of dilations (or erodes), def=1
%  erodeflag - erode instead of dilate
%  flag2d - perform dilation or erosion in 2d only.
%
% Algorithm: a single dilation consists of setting a voxel to 1 if
% it has a neighbor that is 1. The neighborhood is 3x3x3, unless
% flag2d=1, then 3x3x1. Dilation and erosion are reversible only
% if there are no edge effects.
%
% This is pretty slow for big images. Gets slower with the use of
% circshift. The slowness has to do with making a shifted copy.
%


%
% fast_dilate.m
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

md = [];

if(nargin < 1 | nargin > 4)
  fprintf('md = di1(m,<ndilations>,<erodeflag>,<flag2d>)\n');
  return;
end

if(~exist('ndilations','var')) ndilations = []; end
if(isempty(ndilations)) ndilations = 1; end

if(~exist('erodeflag','var')) erodeflag = []; end
if(isempty(erodeflag)) erodeflag = 0; end

if(~exist('flag2d','var')) flag2d = []; end
if(isempty(flag2d)) flag2d = 0; end

if(erodeflag)
  md = ~fast_dilate(~m,ndilations,0,flag2d);
  return;
end

if(ndilations ~= 1)
  m = fast_dilate(m,ndilations-1,erodeflag,flag2d);
end
md = m;

if(flag2d) dsrange = 0;
else       dsrange = [-1:1];
end

md = m;
%tic;
nth = 1;
%fprintf('--------------------\n');
for dc = -1:1
  for dr = -1:1
    for ds = dsrange
      if(dc==0 & dr==0 & ds == 0) continue; end
      vshift = [dc dr];
      if(length(size(m)) > 2)
	vshift = [vshift ds];
      end
      %fprintf('%2d %2d %2d %2d  %g\n',nth,vshift,toc);
      md = md | fast_mshift(m,vshift);
      nth = nth + 1;
    end
  end
end

return;
