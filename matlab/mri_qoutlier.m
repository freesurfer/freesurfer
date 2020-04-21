function qout = mri_qoutlier(y)
% qout = mri_qoutlier(y)
% Quartile distance for outlier detection.
% qout is same size as y (which should be nframes-by-nvox)
% qout value is the distance from the first or third quartile, which
% ever is more divided by the interquartile range. If input value
% is between 1st and 3rd, then qout=0; 
% Quartiles are computed across frame.
% The output can be thresholded to do the outlier detection.
% Typical threshold values are 1.5 or 2
%

if(nargin ~= 1)
  fprintf('qout = mri_qoutlier(y)\n');
  return;
end

nframes = size(y,1);
nvox = size(y,2);

[q1 q2 q3 q4 iqr] = mri_quartiles(y);

q1m = repmat(q1,[nframes 1]);
q3m = repmat(q3,[nframes 1]);
iqrm = repmat(iqr,[nframes 1]);
a = max(max(q1m-y,y-q3m),0);
qout = a./iqrm;

%lfence = repmat(q1-k*iqr,[nframes 1]);
%ufence = repmat(q3+k*iqr,[nframes 1]);
%qout = (y < lfence) | (y > ufence);
%[y q1m q3m iqrm lfence ufence a b qout b>k]
%keyboard

return;

