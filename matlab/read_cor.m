function a = read_cor(fname)
% a = read_cor(fname)


%
% read_cor.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
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



a = zeros(256, 256, 256);

for i=1:256
%	disp(i);
	fname1 = sprintf('%s/COR-%03d', fname, i);
	fid = fopen(fname1, 'r');
	a(:, :, i) = fread(fid, [256 256], 'uchar');
	fclose(fid);
end

% eof
