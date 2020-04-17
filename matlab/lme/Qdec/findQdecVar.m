function col = findQdecVar(Qdec,name)
% col = findQdecVar(Qdec,name)
%
% Gives the column of a variable in a Freesurfer's Qdec table that was read 
% into the two dimensional cell string array Qdec. Returns empty if the 
% variable is not in the Qdec table.
%
% Input
% Qdec: Two dimensional cell string array of Qdec data (eg. read with 
% fReadQdec).
% name: A string with the name of the variable. 
%
% Output
% col: Column of the variable in the cell string array Qdec.
%
% Original Author: Jorge Luis Bernal Rusiel 
%
if nargin < 2
    error('Too few inputs');
end;
col = find(strcmp(getQdecVars(Qdec),name));
