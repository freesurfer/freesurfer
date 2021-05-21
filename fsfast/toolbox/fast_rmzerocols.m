function [M, indnz, indz] = fast_rmzerocols(M0)
% [M indnz indz] = fast_rmzerocols(M0)
%
% Removes columns from M0 that are all zero
% indnz are the indicies of the non-zero columsn
% indz  are the indicies of the zero columsn


%
% fast_rmzerocols.m
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

if(nargin ~= 1)
  fprintf('[M indnz indz] = fast_rmzerocols(M0)\n');
  return;
end

ncols = size(M0,2);

nth = 1;
M = [];
indnz = [];
indz = [];
for c = 1:ncols
  indnzc = find(M0(:,c) ~= 0);
  if(~isempty(indnzc))
    M = [M M0(:,c)];
    indnz = [indnz c];
    nth = nth + 1;
  else
    indz = [indz c];
  end
end

return;
