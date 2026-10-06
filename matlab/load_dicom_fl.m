function [vol, M, dcminfo, mr_parms] = load_dicom_fl(flist)
% [vol, M, dcminfo, mr_parms] = load_dicom_fl(flist)
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
% mr_parms = [tr flipangle te ti]
%
% Does not handle multiple frames correctly yet.
%


%
% load_dicom_fl.m
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


vol=[];
M=[];

if(nargin ~= 1)
  fprintf('[vol, M] = load_dicom_fl(flist)\n');
  return;
end

nfiles = size(flist,1);

tic
fprintf('Loading dicom info for each file \n');
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

ismultiframe = false;
if numel(dcminfo0)==1
  % check whether it is multiframe, or just a dicome sequence with a single slice
  ismultiframe = isfield(dcminfo0, 'SharedFunctionalGroupsSequence') && isfield(dcminfo0, 'PerFrameFunctionalGroupsSequence');
end

if ~ismultiframe
  % Sort by slice location %
  dcminfo = sort_by_sliceloc(dcminfo0);
  
  ipp1 = dcminfo(1).ImagePositionPatient; % position of the first slice
  ipp2 = dcminfo(2).ImagePositionPatient; % position of the second slice
  ippn = dcminfo(nfiles).ImagePositionPatient; % position of the last slice
  ipo1 = dcminfo(1).ImageOrientationPatient; % orientation information
  ps   = dcminfo(1).PixelSpacing; % pixel spacing within the slice
  nslice = nfiles;

else
  % the dicominfo function from matlab produces an ugly structure, with a
  % lot of 'Item_x' fields, at least in MATLAB2022b, the next line cleans
  % the structure a bit.
  dcminfo = clean_dicominfo(dcminfo0); % there's just a single header
  
  ipp1 = dcminfo.PerFrameFunctionalGroupsSequence(1).PlanePositionSequence.ImagePositionPatient;
  ipp2 = dcminfo.PerFrameFunctionalGroupsSequence(2).PlanePositionSequence.ImagePositionPatient;
  ippn = dcminfo.PerFrameFunctionalGroupsSequence(end).PlanePositionSequence.ImagePositionPatient;
  ipo1 = dcminfo.PerFrameFunctionalGroupsSequence(1).PlaneOrientationSequence.ImageOrientationPatient;
  ps   = dcminfo.PerFrameFunctionalGroupsSequence(1).PixelMeasuresSequence.PixelSpacing;
  nslice = dcminfo.NumberOfFrames;
end

% Slice direction cosine %
sdc = ippn - ipp1;
sdc = sdc /sqrt(sum(sdc.^2));

% Distance between slices %
dslice = sqrt(sum((ipp2 - ipp1).^2));

% Matrix of direction cosines %
Mdc = zeros(3,3);
Mdc(:,1) = ipo1(1:3);
Mdc(:,2) = ipo1(4:6);
Mdc(:,3) = sdc;

% Voxel resolution %
delta = [ps; dslice];
D = diag(delta);

% XYZ of first voxel in first slice %
P0 = ipp1;

% Change Siemens to be RAS %
% GE is also LPS - change it to be RAS, too - ebeth %
Manufacturer = dcminfo(1).Manufacturer;
if(strncmpi(Manufacturer,'Siemens',7) | strncmpi(Manufacturer,'ge medical systems',18))
  % Change to RAS
  Mdc(1,:) = -Mdc(1,:); 
  Mdc(2,:) = -Mdc(2,:); 
  P0(1)    = -P0(1);
  P0(2)    = -P0(2);
end

% Compute vox2ras transform %
M = [Mdc*D P0; 0 0 0 1];
if (0&strcmpi(Manufacturer,'Siemens'))
  % Correcting for P0 being at corner of 
  % the first voxel instead of at the center
  M = M*[[eye(3) [0.5 0.5 0]']; 0 0 0 1];  %'
end

% Pre-allocate vol. Note: column and row designations do
% not really mean anything. The "column" is the fastest
% dimension. The "row" is the next fastest, etc.
ndim1 = dcminfo(1).Columns;
ndim2 = dcminfo(1).Rows;
ndim3 = nslice;
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
  
  if ismultiframe
    vol = squeeze(x);
    vol = permute(vol, [2 1 3]);
  else
    % Note: dicomread will transposed the image. This is supposed
    % to help. Anyway, to make the vox2ras transform agree with
    % the volume, the image is transposed back.
    vol(:,:,n) = x'; %'
  end
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

% Lines below correct for the Z-offset in GE machines - ebeth %
% We're told ge machines recenter along superior/inferior axis but
% don't update c_ras - but now c_s should be zero.
if(strcmpi(Manufacturer,'ge medical systems'))
  % Lines below correct for the Z-offset in GE machines
  firstZ = ipp1(3);
  lastXYZ = M*[size(vol)';1]; %'
  % size(imvol) = number of slices in all 3 dirs
  lastZ = lastXYZ(3);
  offsetZ = (lastZ + firstZ)/2.0;
  % Z0 = Z + offsetZ;  [XYZ1]' = M*[CRS1]', need to add to M(3,4)(?)
  M(3,4) = M(3,4) - offsetZ;
end

% if(strcmpi(Manufacturer,'ge medical systems')) 
%   M(3,4) = 0; % Wow - that was actually completely wrong!
% end

% Pull out some info from the header %
if(isfield(dcminfo(1),'FlipAngle')) FlipAngle = pi*dcminfo(1).FlipAngle/180; 
else FlipAngle = 0;
end
if(isfield(dcminfo(1),'EchoTime')) EchoTime = dcminfo(1).EchoTime; 
else EchoTime = 0;
end
if(isfield(dcminfo(1),'RepetitionTime')) RepetitionTime = dcminfo(1).RepetitionTime; 
else RepetitionTime = 0;
end
InversionTime = 0;
mr_parms = [RepetitionTime FlipAngle EchoTime InversionTime];

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

%----------------------------------------------------------%
function dcminfo_out = clean_dicominfo(dcminfo_in)

% subfunction to remove all instances from Item_x as fields from a
% (sub-)structure

fn = fieldnames(dcminfo_in);
if all(strncmp(fn, 'Item_', 5))
  for k = 1:numel(fn)
    dcminfo_out(k) = dcminfo_in.(fn{k});
  end
else
  for k = 1:numel(fn)
    tmp = dcminfo_in.(fn{k});
    if isstruct(tmp) && numel(fieldnames(tmp))>0 && all(strncmp(fieldnames(tmp), 'Item_', 5))
      dcminfo_out.(fn{k}) = clean_dicominfo(tmp);
      
      % recurse
      if isstruct(dcminfo_out.(fn{k}))
        for kk = 1:numel(dcminfo_out.(fn{k}))
          dcminfo_out.(fn{k})(kk) = clean_dicominfo(dcminfo_out.(fn{k})(kk));
        end
      end
    elseif isstruct(tmp) && numel(fieldnames(tmp))>0
      fn_tmp = fieldnames(tmp);
      dcminfo_out.(fn{k}) = struct;
      for kk = 1:numel(fn_tmp)
        tmp2 = tmp.(fn_tmp{kk});
        if isstruct(tmp2)
          dcminfo_out.(fn{k}).(fn_tmp{kk}) = clean_dicominfo(tmp2); 
        else
          dcminfo_out.(fn{k}).(fn_tmp{kk}) = tmp2;
        end
      end
    else
      dcminfo_out.(fn{k}) = tmp;
    end
  end
end

%----------------------------------------------------------%

