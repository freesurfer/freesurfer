function [vertices, label, colortable] = Read_Brain_Annotation(filename)

% [vertices, label, colortable] = Read_Brain_Annotation(annotfilename.annot)
%
% vertices expected to be simply from 0 to number of vertices - 1;
% label is the vector of annotation
%
% colortable is empty struct if not embedded in .annot. Else, it will be
% a struct.
% colortable.numEntries = number of Entries
% colortable.orig_tab = name of original colortable
% colortable.struct_names = list of structure names (e.g. central sulcus and so on)
% colortable.table = n x 5 matrix. 1st column is r, 2nd column is g, 3rd column
% is b, 4th column is flag, 5th column is resultant integer values
% calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0.

fp = fopen(filename, 'r', 'b');
A = fread(fp, 1, 'int');

tmp = fread(fp, 2*A, 'int');
vertices = tmp(1:2:end);
label = tmp(2:2:end);

bool = fread(fp, 1, 'int');
if(isempty(bool)) %means no colortable
   disp('No Colortable found.');
   colortable = struct([]);
   fclose(fp);
   return; 
end

if(bool)
    
    %Read colortable
    numEntries = fread(fp, 1, 'int');    
    colortable.numEntries = numEntries;
    
    len = fread(fp, 1, 'int');
    colortable.orig_tab = fread(fp, len, '*char')';
    colortable.orig_tab = colortable.orig_tab(1:end-1);
    
    colortable.struct_names = cell(numEntries,1);
    colortable.table = zeros(numEntries,5);
    for i = 1:numEntries
       len = fread(fp, 1, 'int');
       colortable.struct_names{i} = fread(fp, len, '*char')';
       colortable.struct_names{i} = colortable.struct_names{i}(1:end-1);
       colortable.table(i,1) = fread(fp, 1, 'int');
       colortable.table(i,2) = fread(fp, 1, 'int');
       colortable.table(i,3) = fread(fp, 1, 'int');
       colortable.table(i,4) = fread(fp, 1, 'int');
       colortable.table(i,5) = colortable.table(i,1) + colortable.table(i,2)*2^8 + colortable.table(i,3)*2^16 + colortable.table(i,4)*2^24;
    end
    disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
    
else
    disp('Error! Should not be expecting bool = 0');    
end

fclose(fp);


