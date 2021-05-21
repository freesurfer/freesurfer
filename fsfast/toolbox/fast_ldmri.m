function [fastvol, msg] = fast_ldmri(volid, volformat, rows, cols, slices, planes)
% fastvol = fast_ldmri(volid, volformat, rows, cols, slices, planes)


%
% fast_ldmri.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

if(~exist('rows'))   rows   = []; end; 
if(~exist('cols'))   cols   = []; end;
if(~exist('slices')) slices = []; end;
if(~exist('planes')) planes = []; end;

[fastvol msg] = fast_ldbvolume(volid, rows, cols, slices, planes);


return;
