function [vol, M, dcminfo] = load_dicom_fl(flist)
% [vol, M, dcminfo] = load_dicom_fl(flist)
%
% Loads a volume from the dicom files in flist.
%
% The volume dimensions are arranged such that the
% readout dimension is first, followed by the phase-encode,
% followed by the slices (this is not implemented yet)
%
% M is the 4x4 vox2ras transform such that
% vol(i1,i2,i3), xyz1 = M*[i1 i2 i3 1] where the
% indicies are 0-based. 
%
% Does not handle multiple frames correctly yet.
%
% $Id: load_dicom_fl.m,v 1.5 2003/01/23 19:57:34 greve Exp $

vol=[];
M=[];

if(nargin ~= 1)
  fprintf('[vol, M] = load_dicom_fl(flist)\n');
  return;
end

nfiles = size(flist,1);

tic
fprintf('Loading dicom info foreach file \n');
for n = 1:nfiles
  fname = deblank(flist(n,:));
  % fprintf('n = %d/%d, %s   %g\n',n,nfiles,fname,toc);
  tmpinfo = dicominfo(fname);
  if(isempty(tmpinfo)) 
    fprintf('ERROR: reading %s\n',fname);
    return;
  end
  if(n > 1)
    % Check that the nth series number agrees with the first
    if(tmpinfo.SeriesNumber ~= dcminfo0(1).SeriesNumber)
      fprintf('ERROR: series number inconsistency (%s)\n',fname);
      return;
    end
  end
  tmpinfo.fname = fname;
  dcminfo0(n) = tmpinfo;
end

% Sort by slice location %
dcminfo = sort_by_sliceloc(dcminfo0);

% Slice direction cosine %
sdc = dcminfo(nfiles).ImagePositionPatient-dcminfo(1).ImagePositionPatient;
sdc = sdc /sqrt(sum(sdc.^2));

% Distance between slices %
dslice = sqrt(sum((dcminfo(2).ImagePositionPatient-dcminfo(1).ImagePositionPatient).^2));

% Matrix of direction cosines %
Mdc = zeros(3,3);
Mdc(:,1) = dcminfo(1).ImageOrientationPatient(1:3);
Mdc(:,2) = dcminfo(1).ImageOrientationPatient(4:6);
Mdc(:,3) = sdc;

% Voxel resolution %
delta = [dcminfo(1).PixelSpacing; dslice];
D = diag(delta);

% XYZ of first voxel in first slice %
P0 = dcminfo(1).ImagePositionPatient;

% Change Siemens to be RAS %
Manufacturer = dcminfo(1).Manufacturer;
if(strcmpi(Manufacturer,'Siemens'))
  % Change to RAS
  Mdc(1,:) = -Mdc(1,:); 
  Mdc(2,:) = -Mdc(2,:); 
  P0(1)    = -P0(1);
  P0(2)    = -P0(2);
  % Correcting for P0 being at corner of 
  % the first voxel instead of at the center
  M = [Mdc*D P0; 0 0 0 1];
  M = M*[[eye(3) [0.5 0.5 0]']; 0 0 0 1];  %'
else
  M = [Mdc*D P0; 0 0 0 1];
end

% Pre-allocate vol. Note: column and row designations do
% not really mean anything. The "column" is the fastest
% dimension. The "row" is the next fastest, etc.
ndim1 = dcminfo(1).Columns;
ndim2 = dcminfo(1).Rows;
ndim3 = nfiles;
vol = zeros(ndim1,ndim2,ndim3);

fprintf('Loading data from each file.\n');
for n = 1:nfiles
  %fprintf('n = %d, %g\n',n,toc);
  fname = dcminfo(n).fname;
  x = dicomread(fname);
  if(isempty(x))
    fprintf('ERROR: could not load pixel data from %s\n',fname);
    return;
  end
  % Note: dicomread will transposed the image. This is supposed
  % to help. Anyway, to make the vox2ras transform agree with
  % the volume, the image is transposed back.
  vol(:,:,n) = x'; %'
end

% Reorder dimensions so that ReadOut dim is first %
if(0 & ~strcmpi(dcminfo(1).PhaseEncodingDirection,'ROW'))
  % This does not work
  fprintf('INFO: permuting vol so that ReadOut is first dim\n');
  vol = permute(vol,[2 1 3]);
  Mtmp = M;
  M(:,1) = Mtmp(:,2);
  M(:,2) = Mtmp(:,1);
end

return;

%----------------------------------------------------------%
function dcminfo2 = sort_by_sliceloc(dcminfo)
  nslices = length(dcminfo);
  sliceloc = zeros(nslices,1);
  for n = 1:nslices
    sliceloc(n) = dcminfo(n).SliceLocation;
  end

  [tmp ind] = sort(sliceloc);
  dcminfo2 = dcminfo(ind);
return;

%----------------------------------------------------------%

