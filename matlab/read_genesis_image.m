function m = read_genesis_image(fname)
% m = read_ge(fname)


%
% read_genesis_image.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
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

fid = fopen(fname, 'r', 'b');
if fid < 0
	disp(sprintf('read_genesis_image: error opening file %s', fname));
	return;
end

fseek(fid, 4, 'bof');
image_offset = fread(fid, 1, 'int');
rows = fread(fid, 1, 'int');
cols = fread(fid, 1, 'int');

fseek(fid, image_offset, 'bof');
m = fread(fid, [rows cols], 'short');

fclose(fid);

% eof
