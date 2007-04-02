function [dat] = load_csd(csdfile)
% dat = load_csd(csdfile)
% 
% reads in Cluster Simulation Data (CSD) as produced by mri_glmfit.
% dat has 5 columes:
%  1. Row number (0-based)
%  2. nClusters 
%  3. MaxClustSize
%  4. MaxSig    
%  5. MaxStat
% 
% Currently does not read in the CSD header.
%
% $Id: load_csd.m,v 1.1 2007/04/02 22:56:02 greve Exp $

%
% load_csd.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/04/02 22:56:02 $
%    $Revision: 1.1 $
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
  fprintf('dat = load_csd(csdfile)\n');
  return;
end

fid = fopen(csdfile);
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
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  dat(nthrow,:) = sscanf(tline,'%f',5);
  nthrow = nthrow + 1;
end % while (1)

fclose(fid);






