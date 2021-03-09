function F = flac_tfilter(flac)
% F = flac_tfilter(flac)
% Builds a matrix to implement a temporal filter. flac.ntp and
% flac.TR must be set. See flac_tfilter_parse.m.
%

%
% flac_tfilter.m
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


  
