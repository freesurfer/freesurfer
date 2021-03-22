function [rg, cg] = subbrain2ghost(imgsize,rb,cb,pedim)
% [rg, cg] = subbrain2ghost(imgsize,rb,cb,<pedim>)
%
% Computes the row and col of a ghost (rg,cg) from the row and
% col of the brain (ie, main image). Note that it can be used
% in the opposite direction (ie, to find the row and col in the
% main image from that of the gost). In fact, running it with
% rg,cg as input will give rb,cb as a result. The ghost of the
% ghost is the orginal image.
%
% imgsize = [nrows ncols]
%
% pedim = 2 means that the phase encode changes from 
% one row to the next. This is the default if pedim is
% unspecfied or null.
%
% pedim = 1 means that the phase encode changes from 
% one col to the next.
%
% Beware of matlab's column-major when assigning PE dim.
%
% See also: indbrain2ghost
%
%


%
% subbrain2ghost.m
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

rg = [];
cg = [];

if(nargin < 3 | nargin > 4)
  fprintf('[rg, cg] = subbrain2ghost(imgsize,rb,cb,<pedim>)\n');
  return;
end

if(~exist('pedim','var')) pedim = []; end
if(isempty(pedim)) pedim = 2; end
if(pedim ~= 1 & pedim ~= 2)
  fprintf('ERORR: pedim = %d, must be 1 or 2\n',pedim);
  return;
end

nr = imgsize(1);
nc = imgsize(2);

% Make sure none of the rows, cols, or slices are out of range.
if(length(find(rb>imgsize(1))))
  fprintf('ERROR: rb values are out of range\n');
  return;
end
if(length(find(cb>imgsize(2))))
  fprintf('ERROR: cb values are out of range\n');
  return;
end

if(pedim == 2)
  % phase encode changes from one row to the next.
  rg = rb+round(nr/2);
  indtmp = find(rg>nr);
  rg(indtmp) = rg(indtmp) - nr;
  cg = cb;
  return;
end

if(pedim == 1)
  % phase encode changes from one col to the next.
  cg = cb+round(nc/2);
  indtmp = find(cg>nc);
  cg(indtmp) = cg(indtmp) - nc;
  rg = rb;
  return;
end

% Should never get here
return;


