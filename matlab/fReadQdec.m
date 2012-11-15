function Qdec = fReadQdec(fname)
% Qdec = fReadQdec(fname) 
%
% Reads a Freesurfer's Qdec table file into cell string array Qdec.
%
% Input
% fname: The name of a Qdec table data file. 
%
% Output
% Qdec: Two dimensional cell string array with the data.
%
% $Revision: 1.1 $  $Date: 2012/11/15 15:17:51 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: vinke $
%    $Date: 2012/11/15 15:17:51 $
%    $Revision: 1.1 $
%
fid = fopen(fname);
tline = fgetl(fid);
Qdec = [];
i = 1;
while ischar(tline)
    j = 1;
    [str,remain] = strtok(tline, ' ');   
    while ~isempty(str)
        Qdec{i,j} = str;
        j = j + 1;
        [str,remain] = strtok(remain, ' ');
    end
    i = i + 1;
    tline = fgetl(fid);
end
fclose(fid); 
