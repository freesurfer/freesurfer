function par = fmri_ldpar(varargin)
%
% Load specified par files for multiple runs.
%
% par = fmri_ldpar(ParFileList)
% par = fmri_ldpar(ParFile1, ParFile2, ...)
%
% ParFileList is a vertical cat of parfile names (ie, each run's'
% parfile name on a different row).
%
% ParFile format is as follows:
%   Column 1: stimulus presentation time (float or int)
%   Column 2: stimulus condition number
%   Column > 2: ignored
%  
% par dimensionality: nPresentations x 2 x nRuns
%
% $Id: fmri_ldpar.m,v 1.1 2003/03/04 20:47:39 greve Exp $

par = [];

if(nargin == 0)
  msg = 'USAGE: par = fmri_ldpar(ParFile)';
  qoe(msg);
  error(msg);
end

if( length(varargin) == 1)
  ParFile = varargin{1};
  nRuns = size(ParFile,1);
else
  nRuns = length(varargin);
  ParFile = '';
  for r = 1:nRuns,
    ParFile = strvcat(ParFile,varargin{r});
  end
end

% Go through each run and load par %
for r=1:nRuns,

  %%% Open the par file %%%%
  [fid msg] = fopen(deblank(ParFile(r,:)),'r');
  if fid == -1 
    fprintf('%s\n',msg);
    fprintf('%s\n',pwd);
    msg = sprintf('Could not open %s',ParFile(r,:)); 
    qoe(msg) ; error(msg);
  end

  % Read all the lines %
  sLines = readln(fid);
  fclose(fid);

  if(isempty(sLines)) 
    fprintf('ERROR: reading parfile %s\n',ParFile(r,:));
    return; 
  end

  % Go through each line %
  for n = 1:size(sLines,1)
    tmp = sscanf(sLines(n,:),'%f',2)'; %'
    if(~isempty(tmp))
      par(n,:,r) = tmp;
    end
    % par(n,:,r) = sscanf(sLines(n,:),'%f',2)'; %'
  end

end

return;


