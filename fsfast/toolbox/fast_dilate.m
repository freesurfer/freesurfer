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
% $Id: fast_dilate.m,v 1.1 2005/02/17 18:19:16 greve Exp $

md = [];

if(nargin < 1 | nargin > 4)
  fprintf('md = fast_dilate(m,<ndilations>,<erodeflag>,<flag2d>)\n');
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
else       dsrange = 1;
end


for r = 1:size(m,1)
  for c = 1:size(m,2)
    for s = 1:size(m,3)
      if(md(r,c,s)) continue; end

      for dr = -1:1
	rn = r+dr;
	if(rn < 1 | rn > size(m,1)) continue; end
	for dc = -1:1
	  cn = c+dc;
	  if(cn < 1 | cn > size(m,2)) continue; end
	  for ds = -dsrange:dsrange
	    sn = s+ds;
	    if(sn < 1 | sn > size(m,3)) continue; end
	    if(m(rn,cn,sn))
	      md(r,c,s) = 1;
	      continue;
	    end
	  end % ds
	end % dc
      end % dr
    
    end
  end
end

return;

