function synch = fmri_ldsynch(synchfile)
%
% synch = fmri_ldsynch(synchfile)
%
% Load specified slice synchronization file
%
% $Id: fmri_ldsynch.m,v 1.1 2003/03/04 20:47:40 greve Exp $

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


