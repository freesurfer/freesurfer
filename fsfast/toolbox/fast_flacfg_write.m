function err = fast_flacfg_write(c,cfgfile)
% err = fast_flacfg_write(cfg,cfgfile)


%
% fast_flacfg_write.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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


err = 1;

if(nargin ~= 2)
  fprintf('err = fast_flacfg_write(cfg,cfgfile)\n');
  return;
end

fid = fopen(cfgfile,'w');
if(fid == -1)
  fprintf('ERROR: cannot open cfg file %s\n',cfgfile);
  return;
end

fprintf(fid,'FSFAST-FLACFG %d\n',c.version);
if(~isempty(c.flaname)) fprintf(fid,'flaname %s\n',c.flaname); end
if(~isempty(c.TR)) fprintf(fid,'tr %f\n',c.TR); end
if(~isempty(c.inorm)) fprintf(fid,'inorm %f\n',c.inorm); end
if(~isempty(c.nskip)) fprintf(fid,'nskip %d\n',c.nskip); end
if(~isempty(c.slicetiming)) fprintf(fid,'slicetiming %s\n',c.slicetiming); end
if(~isempty(c.evschfname)) fprintf(fid,'evschfname %s\n',c.evschfname); end
if(~isempty(c.runlistfile)) fprintf(fid,'runlistfile %s\n',c.runlistfile); end
if(~isempty(c.fsd)) fprintf(fid,'fsd %s\n',c.fsd); end
if(~isempty(c.funcstem)) fprintf(fid,'funcstem %s\n',c.funcstem); end
if(~isempty(c.maskstem)) fprintf(fid,'maskstem %s\n',c.maskstem); end
if(~isempty(c.usetpexclude)) fprintf(fid,'usetpexclude %d\n',c.usetpexclude); end

nfx = length(c.fxlist);
for nthfx = 1:nfx
  c.nthfx = nthfx;
  tline = fast_fxcfg('createline',c);
  fprintf(fid,'%s\n',tline);
end

fclose(fid);

err = 0;
return;


