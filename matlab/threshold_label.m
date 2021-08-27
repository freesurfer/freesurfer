function err = threshold_label(sname, lname, thresh, outname)
% err = threshold_label(sname, lname, thresh, outname)
%
% reads the label file 'lname' (do not include .label extension) from
% the subject 'sname' in the subject's label directory, removes any
% nodes whose stat is less than thresh, saves in outname.
% 
% Example:
% threshold_label('mysubject','lh.V1',.8,'lh.V1.thresh.label');

%
% threshold_label.m
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


err = 0;
if(nargin ~= 4)
  fprintf('threshold_label(sname, lname, thresh, outname)\n');
  return;
end

l = read_label(sname,lname);
if(isempty(l)) err=1; return; end

stat = l(:,5);
ind = find(stat > thresh);
lindex = l(ind,1);
lxyz = l(ind,[2:4]);
lvals = l(ind,5);
ok = write_label(lindex, lxyz, lvals, outname, sname);

return;

