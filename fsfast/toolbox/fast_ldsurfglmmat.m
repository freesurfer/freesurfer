function [subjnames, desmat] = fast_ldsurfglmmat(fname)
% [subjnames, desmat] = fast_ldsurfglmmat(fname)


%
% fast_ldsurfglmmat.m
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

subjnames = [];
desmat = [];

if(nargin ~= 1)
  fprintf('USAGE: [subjnames, desmat] = fast_ldsurfglmmat(fname)\n');
  return;
end

fid = fopen(fname);
if(fid == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

% Count the number of items %
nitems = 0;
while(1)
  s = fscanf(fid,'%s',1);
  if(isempty(s)) break; end
  nitems = nitems + 1;
  %fprintf('%d %s\n',nitems,s);
end
fprintf('INFO: Found %d items in matrix\n',nitems);

% Close and open to return to start
fclose(fid);
fid = fopen(fname);

% Count the number of lines %
nlines = 0;
while 1
  tline = fgetl(fid);
  if ~ischar(tline) break; end
  nlines = nlines + 1;
end
fprintf('INFO: Found %d lines in matrix\n',nlines);

% Close and open to return to start
fclose(fid);
fid = fopen(fname);

nsubjects = nlines;
ncols = nitems/nlines - 1;

% Read in the data
fmt = ['%*s ' repmat('%f ',[1 ncols])];
nlines = 0;
while 1
  tline = fgetl(fid);
  if ~ischar(tline) break; end
  nlines = nlines + 1;
  subj = sscanf(tline,'%s',1);
  subjnames = strvcat(subjnames,subj);
  d = sscanf(tline,fmt,ncols)'; %'
  desmat = [desmat; d];
end

fclose(fid);

return;
