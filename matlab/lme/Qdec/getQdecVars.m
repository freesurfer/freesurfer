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
% Original Author: Jorge Luis Bernal Rusiel 
%
Vars = Qdec(1,:)';
