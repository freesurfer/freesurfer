function [dat thresh fwhm searchspace anattype] = load_csd(csdfile)
% [dat thresh fwhm searchspace anattype] = load_csd(csdfile)
% 
% reads in Cluster Simulation Data (CSD) as produced by mri_glmfit.
% dat has 5 columes:
%  1. Row number (0-based)
%  2. nClusters 
%  3. MaxClustSize
%  4. MaxSig    
%  5. MaxStat
% 
% thresh is -log10(p)
% fwhm in mm^D
% searchspace in mm^D
% anattype is volume or surface
%
% $Id: load_csd.m,v 1.4 2011/03/02 00:04:12 nicks Exp $

%
% load_csd.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
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

if(nargin ~= 1)
  nargin
  fprintf('dat = load_csd(csdfile)\n');
  return;
end


%'synth-glm-surf/n0001/csd/mc-z.pos.j001-osgm.csd'
fid = fopen(csdfile,'r');
if(fid == -1)
  fprintf('ERROR: opening %s\n',csdfile);
  return;
end

tline = fgetl(fid);
if(tline == -1)
  fprintf('ERROR: %s is not correctly formatted, no first line\n', ...
	  csdfile);
  fclose(fid);
  return;
end

%----------- Loop through all the lines ----------------------%
nthrow = 1;
while(1)

  % scroll through any blank lines or comments %
  while(1)
    tline = fgetl(fid);
    if(~isempty(tline))
      if(tline(1) ~= '#') break; end
      key = sscanf(tline,'%*s %s',1);
      if(strcmp(key,'thresh'))
	thresh = sscanf(tline,'%*s %*s %f');
      end
      if(strcmp(key,'nullfwhm'))
	fwhm = sscanf(tline,'%*s %*s %f');
      end
      if(strcmp(key,'searchspace'))
	searchspace = sscanf(tline,'%*s %*s %f');
      end
      if(strcmp(key,'anattype'))
	anattype = sscanf(tline,'%*s %*s %s');
      end
    end
  end
  if(tline(1) == -1) break; end
  dat(nthrow,:) = sscanf(tline,'%f',5);

  nthrow = nthrow + 1;
end % while (1)

fclose(fid);

return;








