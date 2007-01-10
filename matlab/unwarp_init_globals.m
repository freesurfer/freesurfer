function unwarp_init_globals(called_by_script)


% Called by convert_unwarp_resample.m and unwarp_scanners_table.m
% This can't be the right way to do this.


%
% unwarp_init_globals.m
%
% Original Author: Elizabeth Haley
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision: 1.3 $
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
