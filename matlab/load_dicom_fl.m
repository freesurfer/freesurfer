function [vol, M, dcmvolinfo] = load_dicom_fl(flist)
% [vol, M, dcminfo] = load_dicom_fl(flist)
%
% Loads a volume from the dicom files in flist.

vol=[];
M=[];
dcmvolinfo = [];

if(nargin ~= 1)
  fprintf('[vol, M, dcminfo] = load_dicom_fl(flist)\n');
  return;
end

nfiles = size(flist,1);

tic
for n = 1:nfiles
  fprintf('n = %d, %g\n',n,toc);
  fname = deblank(flist(n,:));
  tmpinfo = dicominfo(fname);
  if(isempty(tmpinfo)) 
    fprintf('ERROR: reading %s\n',fname);
    return;
  end
  if(n > 1)
    if(tmpinfo.SeriesNumber ~= dcminfo(1).SeriesNumber)
      fprintf('ERROR: series number inconsistency (%s)\n',fname);
      return;
    end
  end
  tmpinfo.fname = fname;
  dcminfo(n) = tmpinfo;
end

dcminfo2 = sort_by_sliceloc(dcminfo);

sdc = dcminfo(nfiles).ImagePositionPatient-dcminfo(1).ImagePositionPatient;
sdc = sdc /sqrt(sum(sdc.^2));

keyboard
return;

%----------------------------------------------------------%
function dcminfo2 = sort_by_sliceloc(dcminfo)
  nslices = length(dcminfo);
  sliceloc = zeros(nslices,1);
  for n = 1:nslices
    sliceloc = dcminfo(n).SliceLocation;
  end

  [tmp ind] = sort(sliceloc);
  dcminfo2 = dcminfo(ind);
return;

%----------------------------------------------------------%



%----------------------------------------------------------%
