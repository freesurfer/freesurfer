function sfa = fmri_ldsfa(stem)
% sfa = fmri_ldsfa(stem)


%
% fmri_ldsfa.m
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

if(nargin ~= 1) 
  msg = 'USAGE: fmri_ldsfa(stem)';
  qoe(msg); error(msg);
end

stem = deblank(stem);

%% Open the sfa file %%
sfaname = sprintf('%s.sfa',stem);
fid = fopen(sfaname,'r');
if(fid == -1) 
  msg = sprintf('Could not open %s for reading',sfaname);
  qoe(msg); error(msg);
end

%% Get past the first line %%
dummy = fscanf(fid,'%s',1);

sfa = ldstruct(fid);
fclose(fid);
return;



%% Read the rest of the lines %%
slines = readln(fid);

%% Set up the structure %%
sfa = fmri_sfastruct;

%% Go through each line and get the value for the strcture %%
for n = 1:size(slines,1)
  fld = sscanf(slines(n,:),'%s',1)

  switch(fld)
    case 'isignal',
      fsa.isgnal = fscanf(fid,'%d',fsa.nsignal);

    case 'inoise',
      fsa.inoise = fscanf(fid,'%d',fsa.nnoise);

    case 'iexclude',
      fsa.iexclude = fscanf(fid,'%d',fsa.nexclude);

    case 'slice_delay',
      fsa.slice_delay = fscanf(fid,'%f',fsa.nslices);

    otherwise,
       

  end %- switch -%


end



fclose(fid);
return;
