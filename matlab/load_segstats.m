function [segname, segindex, segstats] = load_segstats(segstatsfile,subject)
% [segname, segindex, segstats] = load_segstats(segstatsfile,<subject>)
%
% Loads data from stats file produced by mri_segstats. This
% includes the aseg.stats produced by recon-all.
%
% segstatsfile - path to seg stats file
% subject - optional. Looks for segstatsfile relative to 
%    SUBJECTS_DIR/subject/stats. If segstatsfile is empty,
%    then loads SUBJECTS_DIR/subject/stats/aseg.stats
%
%  segname is a list of the segmentation names 
%    (col 5 from aseg.stats file)
%  segindex is the index of each segmentation
%  segstats are the stats for each index. Each col is a stat:
%    1. number of voxels
%    2. volume of voxels (mm^3) -- same as number but scaled by voxvol
%    3. mean intensity over space
%    4. std intensity over space
%    5. min intensity over space
%    6. max intensity over space
%    7. range intensity over space
%
% The intensities are from whatever was passed as the --in volume
% to mri_segstats. For aseg.stats, this is the norm.mgz volume.
%
% Any line that begins with a '#' is ignored.
%
% Examples:
% [segname segindex segstats] = load_segstats(segstatsfile);
% [segname segindex segstats] = load_segstats(segstatsfile,'bert');
% [segname segindex segstats] = load_segstats([],'bert');
%
%


%
% load_segstats.m
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

segname='';
segindex=[];
segstats=[];

if(nargin < 1 |  nargin > 2)
  fprintf('[segname, segindex, segstats] = load_segstats(segstatsfile,<subject>)\n');
  return;
end

if(~exist('subject','var')) subject = []; end
if(~isempty(subject))
  if(isempty(segstatsfile)) segstatsfile = 'aseg.stats'; end
  SUBJECTS_DIR = getenv('SUBJECTS_DIR');
  segstatsfile = sprintf('%s/%s/stats/%s',SUBJECTS_DIR,subject,segstatsfile);
end

fid = fopen(segstatsfile);
if(fid == -1)
  fprintf('ERROR: opening %s\n',segstatsfile);
  return;
end

tline = fgetl(fid);
if(tline == -1)
  fprintf('ERROR: %s is not correctly formatted, no first line\n', ...
	  segstatsfile);
  fclose(fid);
  return;
end

%----------- Loop through all the lines ----------------------%
nthrow = 1;
while(1)

  % scroll through any blank lines or comments (#)
  while(1)
    tline = fgetl(fid);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  indx   = sscanf(tline,'%d',1);
  segid  = sscanf(tline,'%*d %d',1);
  nvox   = sscanf(tline,'%*d %*d %d',1);
  vol    = sscanf(tline,'%*d %*d %*d %f',1);
  %segnm  = sscanf(tline,'%*d %*d %*d %*f %s',1);
  segnm  = sscanfitem(tline,5);
  segmn  = sscanf(tline,'%*d %*d %*d %*f %*s %f',1);
  segstd = sscanf(tline,'%*d %*d %*d %*f %*s %*f %f',1);
  segmin = sscanf(tline,'%*d %*d %*d %*f %*s %*f %*f %f',1);
  segmax = sscanf(tline,'%*d %*d %*d %*f %*s %*f %*f %*f %f',1);
  segrng = sscanf(tline,'%*d %*d %*d %*f %*s %*f %*f %*f %*f %f',1);

  segname = strvcat(segname,segnm);
  segindex = [segindex; segid];
  segstats(nthrow,:) = [nvox vol segmn segstd segmin segmax segrng];
  nthrow = nthrow + 1;
end % while (1)

fclose(fid);

return;
%---------------------------------------------------------%
