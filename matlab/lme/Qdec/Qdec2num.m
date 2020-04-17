function M = Qdec2num(Qdec)
% M = Qdec2num(Qdec)
%
% Attempts to convert the two dimensional cell string array Qdec to numeric
% matrix M.
%
% Input
% Qdec: Two dimensional cell string array of Qdec data (eg. read with 
% fReadQdec).
%
% Output
% M: Numeric matrix (without the name of the variables).
%
% Original Author: Jorge Luis Bernal Rusiel 
%
Qdec = Qdec(2:end,:);
M = str2double(Qdec);
