function a = read_cor(fname)
% a = read_cor(fname)


%
% read_cor.m
%
% Original Author: Bruce Fischl
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



a = zeros(256, 256, 256);

for i=1:256
%	disp(i);
	fname1 = sprintf('%s/COR-%03d', fname, i);
	fid = fopen(fname1, 'r');
	a(:, :, i) = fread(fid, [256 256], 'uchar');
	fclose(fid);
end

% eof
