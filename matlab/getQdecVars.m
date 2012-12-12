function Vars = getQdecVars(Qdec)
% Vars = getQdecVars(Qdec)
%
% Returns a one dimensional cell string array with the name of the variables
% in Qdec.
%
% Input
% Qdec: Two dimensional cell string array of Qdec data.
%
% Output
% Vars: One dimensional cell string array.
%
% $Revision: 1.2 $  $Date: 2012/12/12 22:58:12 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: vinke $
%    $Date: 2012/12/12 22:58:12 $
%    $Revision: 1.2 $
%
Vars = Qdec(1,:)';