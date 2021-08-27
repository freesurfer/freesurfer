function [RescaleFactor, MeanVal] = fast_rescalefactor(MeanValFile, RescaleTarget)
% fast_rescalefactor(MeanValFile, RescaleTarget)


%
% fast_rescalefactor.m
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

if(nargin ~= 2)
  msg = 'USAGE: [factor meanval] = fast_rescalefactor(MeanValFile, RescaleTarget)';
  qoe(msg);error(msg);
end

fid = fopen(MeanValFile);
if(fid == -1)
  msg = sprintf('ERROR: could not open %s',MeanValFile);
  qoe(msg);error(msg);
end

MeanVal = fscanf(fid,'%f');
RescaleFactor = RescaleTarget/MeanVal;

fclose(fid);

return;






