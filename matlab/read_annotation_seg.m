function [seg, ctab] = read_annotation_seg(annotname)
% [seg, ctab] = read_annotation_seg(annotname)
%
% read the annotation as a "seg", ie, a row vector where the value is the
% index into the ctab. This is easier to work with than the silly annot
% format.
%

[vertices label ctab] = read_annotation(annotname);

seg = zeros(1,length(vertices));
for n = 1:ctab.numEntries
  ind = find(label == ctab.table(n,5));
  seg(1,ind) = n-1;
end

% stgctab = strmatch('superiortemporal',char(ctab.struct_names));
% stgcode = ctab.table(stgctab,5);
% indstg = find(label==stgcode);
% nstg = length(indstg);

return



