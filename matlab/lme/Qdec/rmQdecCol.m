function Qdec2 = rmQdecCol(Qdec1,col)
% Qdec2 = rmQdecCol(Qdec1,col)
%
% Removes the specify column from cell string array Qdec1. 
%
% Input
% Qdec1: Two dimensional cell string array of Qdec data (eg. read with 
% fReadQdec).
% col: Column to be removed.
%
% Output
% Qdec2: Two dimensional cell string array of Qdec data.
%
% $Revision: 1.2 $  $Date: 2015/01/06 17:14:53 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: mreuter $
%    $Date: 2015/01/06 17:14:53 $
%    $Revision: 1.2 $
%
if nargin < 2
    error('Too few inputs');
end;
Qdec2 = [Qdec1(:,1:col-1) Qdec1(:,col+1:end)]; 


