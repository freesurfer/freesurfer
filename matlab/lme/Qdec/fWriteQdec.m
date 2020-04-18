function fWriteQdec(fname,Qdec)
% fWriteQdec(fname,Qdec)
%
% Writes cell string array Qdec to a Freesurfer's Qdec table data file.
%
% Input
% fname: The name of a Qdec table data file. 
% Qdec: Two dimensional cell string array with the data.
%
% Original Author: Jorge Luis Bernal Rusiel 
%
if nargin < 2
    error('Too few inputs');
end;
szQ = size(Qdec);
fid = fopen (fname,'w');
if fid < 0
    error('Do you have write permissions for %s?', pwd);
end
if ~strcmp(Qdec{1,1},'fsid')
      warning('Dat{1,1} is not ''fsid'' will change it ');
      Qdec{1,1} = 'fsid';
end;
for i=1:szQ(1)
    j = 1;
    while j < szQ(2)
        fprintf(fid, [char(Qdec{i,j}) ' ']);
        j = j + 1;
    end;
    fprintf(fid, [char(Qdec{i,j}) '\n']);
end;
fclose (fid);
