function unwarp_init_globals(called_by_script)


% Called by convert_unwarp_resample.m and unwarp_scanners_table.m
% This can't be the right way to do this.


%
% unwarp_init_globals.m
%
% Original Author: Elizabeth Haley
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.4 $
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



global GRADWARPPATH;
global TABLE;
global QuitOnError;

GRADWARPPATH = sprintf('%s/grad_unwarp_tables',getenv('DEV'));
TABLE = sprintf('%s/table.mat',GRADWARPPATH);
if (exist('called_by_script') & ~isempty(called_by_script))
  QuitOnError = called_by_script;
else
  QuitOnError = 1;
end

return;
