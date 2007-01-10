function [fastvol, msg] = fast_ldmri(volid, volformat, rows, cols, slices, planes)
% fastvol = fast_ldmri(volid, volformat, rows, cols, slices, planes)


%
% fast_ldmri.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

if(~exist('rows'))   rows   = []; end; 
if(~exist('cols'))   cols   = []; end;
if(~exist('slices')) slices = []; end;
if(~exist('planes')) planes = []; end;

[fastvol msg] = fast_ldbvolume(volid, rows, cols, slices, planes);


return;
