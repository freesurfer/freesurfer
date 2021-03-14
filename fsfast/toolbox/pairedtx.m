function [X, C] = pairedtx(npairs)
% [X C] = pairedtx(npairs)
%
% Creates a design matrix for a paired test for npairs of data. It
% is assumed that each pair appears consecutively. The design
% matrix will have size 2*npairs by npairs+1. C can be used to test
% whether the paired difference is equal to zero. C is a
% vector of length npairs+1 with all the components equaling 0
% except the last which is 1.
%
%


%
% pairedtx.m
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
  fprintf('[X C] = pairedtx(npairs)\n');
  return;
end

S = repmat([1 1]',[1 1 npairs]);
S = fast_blockdiag(S);
D = repmat([0.5 -0.5]',[npairs 1]);
X = [S D];

C = zeros(1,npairs+1);
C(end) = 1;

return;

