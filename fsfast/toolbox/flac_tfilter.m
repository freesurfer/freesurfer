function F = flac_tfilter(flac)
% F = flac_tfilter(flac)
% Builds a matrix to implement a temporal filter. flac.ntp and
% flac.TR must be set. See flac_tfilter_parse.m.
%
% $Id: flac_tfilter.m,v 1.1 2010/12/02 19:13:09 greve Exp $

%
% flac_tfilter.m
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



F = [];
if(nargin ~= 1)
  fprintf('F = flac_tfilter(flac)\n');
  return;
end

tfilter = flac.tfilter;
if(isempty(flac.tfilter)) return; end
  
switch (tfilter.model)
 
 %--------------------------------------------
 case {'lpf'} % lowpass filter
  % 2 parameters: cutoffHz order
  cutoffHz = tfilter.params(1);
  order    = tfilter.params(2);
  F = fast_lpfmtx(cutoffHz,flac.TR,flac.ntp,order);
  
  %--------------------------------------------  
 otherwise
  fprintf('ERROR: flac_tfilter_parse: model %s unrecoginized\n',tfilter.model);
  ev = [];
  return;
  
end % switch
  
  
return;


  
