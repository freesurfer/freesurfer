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


