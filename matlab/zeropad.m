function z = zeropad(a,M,N,val,row0,col0)
%ZEROPAD	Introduces zero-padding in a matrix (or any other value)
%
%	z = ZEROPAD(MAT, M, N, row, col[, val]) pads MAT to a M-by-N matrix
%
%	starting from (row,col) with a value 'val' (default zero)
%
% if only 3 parameters are passed A will be centered

[ma,na] = size(a);
if (nargin < 5)
  row0 = round((M-size(a,1))/2)+1 ;
end
if (nargin < 6)
  col0 = round((N-size(a,2))/2)+1 ;
end

if (nargin < 4)
  val = 0 ;
end
z = val*ones(M,N);

z(row0:row0+ma-1,col0:col0+na-1) = a;

