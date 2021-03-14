function r = fmri_touch(fname);
% r = fmri_touch(fname);
% 
% Simple function to create a file called fname.  This is supposed
% to be something like the unix touch, but it has
% no effect if fname already exists.
%
%


%
% fmri_touch.m
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

if(nargin ~= 1) 
  msg = 'USAGE: r = fmri_touch(fname);';
  qoe(msg); error(msg);
end

fname = deblank(fname);

fid = fopen(fname,'a');
if(fid == -1) 
  msg = sprintf('Could not open %s for appending',fname);
  qoe(msg); error(msg);
end

fclose(fid);

r = 0;
return;
