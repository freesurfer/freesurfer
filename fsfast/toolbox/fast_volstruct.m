function fv = fast_mrivolstruct
% fastvol = fast_mrivolstruct


%
% fast_volstruct.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:05 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

fv.data      = [];
fv.szvol     = []; % not necesarily same as size(fv.data)
fv.szpix     = [];
fv.nvox      = 0;
fv.rows      = [];
fv.cols      = [];
fv.slices    = [];
fv.planes    = [];
fv.volid     = '';
fv.precision = '';
fv.endian    = '';
fv.format    = '';
