function synch = fmri_ldsynch(synchfile)
%
% synch = fmri_ldsynch(synchfile)
%
% Load specified slice synchronization file
%
%


%
% fmri_ldsynch.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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

if(nargin == 0)
  msg = 'USAGE: synch = fmri_ldsynch(SynchFile)';
  qoe(msg);
  error(msg);
end


  %%% Open the par file %%%%
  [fid msg] = fopen(deblank(synchfile),'r');
  if fid == -1 
    fprintf('%s\n',msg);
    fprintf('%s\n',pwd);
    dir(SynchFile)
    qoe(msg);
    error( sprintf('Could not open %s',synchfile)); 
  end

  % Read all the lines %
  sLines = readln(fid);
  fclose(fid);

  % Go through each line %
  for n = 1:size(sLines,1)
    synch(n,:) = sscanf(sLines(n,:),'%f',2)'; %'
  end

return;


