% $Id: unwarp_init_globals.m,v 1.2 2003/07/31 21:47:40 ebeth Exp $
% Called by convert_unwarp_resample.m and unwarp_scanners_table.m
% This can't be the right way to do this.

function unwarp_init_globals(called_by_script)

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
