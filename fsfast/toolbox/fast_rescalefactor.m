function [RescaleFactor, MeanVal] = fast_rescalefactor(MeanValFile, RescaleTarget)
% fast_rescalefactor(MeanValFile, RescaleTarget)


%
% fast_rescalefactor.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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






