function [fastvol, msg] = fast_ldmri(volid, volformat, rows, cols, slices, planes)
% fastvol = fast_ldmri(volid, volformat, rows, cols, slices, planes)

if(~exist('rows'))   rows   = []; end; 
if(~exist('cols'))   cols   = []; end;
if(~exist('slices')) slices = []; end;
if(~exist('planes')) planes = []; end;

[fastvol msg] = fast_ldbvolume(volid, rows, cols, slices, planes);


return;
