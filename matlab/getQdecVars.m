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
% $Revision: 1.2.2.2 $  $Date: 2013/01/22 20:59:08 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/01/22 20:59:08 $
%    $Revision: 1.2.2.2 $
%
Vars = Qdec(1,:)';