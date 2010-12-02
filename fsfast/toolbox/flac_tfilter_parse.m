function tfilter = flac_tfilter_parse(tline)
% tfilter = flac_tfilter_parse(tline)
%
% Parses a TFILTER line from a FLAC file. If a tline is not provided,
% then returns an empty tfilter struct.  When adding filters to this
% file, also add to flac_tfilter.
%
% TFILTER TFilterName <parameters>
%

%
% flac_tfilter_parse.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2010/12/02 19:13:09 $
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

tfilter = [];
if(nargin > 1)
  fprintf('tfilter = fast_tfilter_parse(<tline>)\n');
  return;
end

tfilter.model      = [];  % Name of the model
tfilter.params     = [];  % params
if(nargin == 0)  return; end

% Read in the model type
[tfilter.model c] = sscanf(tline,'%*s %s',1);
if(c ~= 1) fprintf('Format error\n'); tfilter=[]; return; end

switch (tfilter.model)
 
 %--------------------------------------------
 case {'lpf'} % lowpass filter
  % 2 parameters: cutoffHz order
  % TFILTER lpf 0.1 5
  [item c] = sscanfitem(tline,3);
  if(c ~= 1) fprintf('Format error\n'); tfilter=[]; return; end
  tfilter.params(1) = sscanf(item,'%f',1); % cutoffHz
  [item c] = sscanfitem(tline,4);
  if(c ~= 1) fprintf('Format error\n'); tfilter=[]; return; end
  tfilter.params(2) = sscanf(item,'%d',1); % order

  %--------------------------------------------  
 otherwise
  fprintf('ERROR: flac_tfilter_parse: model %s unrecoginized\n',tfilter.model);
  ev = [];
  return;
  
end % switch

return;
