function Qdec = num2Qdec(M,Vars)
% Qdec = num2Qdec(M,Vars)
%
% Converts numeric matrix M to two dimensional cell string array Qdec.
%
% Input
% M: Numeric matrix.
% Vars: One dimensional cell string array with the name of the Qdec 
% variables corresponding to the columns of M.
%
% Output
% Qdec: Two dimensional cell string array of Qdec data .
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:10 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:10 $
%    $Revision: 1.1.2.2 $
%
if nargin < 2
    error('Too few inputs');
end;
szM = size(M);
Qdec = cell(szM(1)+1,szM(2));
Qdec(1,1:end) = Vars;
for i=2:szM(1)+1
    for j=1:szM(2)
        Qdec{i,j} = num2str(M(i-1,j));
    end;
end;