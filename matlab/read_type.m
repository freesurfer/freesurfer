function m = read_type(fname, offset, size, type)
% m = read_type('fname', offset, size, 'type')


%
% read_type.m
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
	disp(sprintf('read_type: error opening file %s', fname));
	return;
end

fseek(fid, offset, 'bof');
m = fread(fid, size, type);

fclose(fid);

% eof
