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
% $Id: pairedtx.m,v 1.1 2005/01/03 23:19:33 greve Exp $

if(nargin ~= 1)
  fprintf('[X C] = pairedtx(npairs)\n');
  return;
end

S = repmat([1 1]',[1 1 npairs]);
S = fast_blockdiag(S);
D = repmat([1 -1]',[npairs 1]);
X = [S D];

C = zeros(1,npairs+1);
C(end) = 1;

return;

