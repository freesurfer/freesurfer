function T = tdr_gradmtx(mask)
% T = tdr_gradmtx(mask)
%
% Computes a matrix to compute a spatial difference between a voxel
% and its 4 nearest neighbors within the mask in a slice.
%
%


%
% tdr_gradmtx.m
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

T = [];

if(nargin ~= 1)
  fprintf('T = tdr_gradmtx(mask)\n');
  return;
end

% Slice indices of voxels in the mask
indm = find(mask);
nm = length(indm);
if(nm == 0)
  fprintf('ERROR: no voxels found in the mask\n');
  return
end
%fprintf('%d voxels found in the mask\n',nm);

% Size of the slice
szm = size(mask);
nmrows = szm(1); % Number of rows in the slice
nmcols = szm(2); % Number of cols in the slice

% Slice rows and cols of the voxels in the mask
[rm cm] = ind2sub(szm,indm);

% Create a slice in which the voxel value is the
% index into the mask voxels.
maskindex = zeros(szm);
maskindex(indm) = [1:nm];

T = zeros(nm,nm);
szT = size(T);

% Go thru each neighbor offset
for dr = [-1 0 1]
  for dc = [-1 0 1]

    % Dont process the center or diagonal neighbors
    if(abs(dr)+abs(dc) > 1) continue; end
    if(dr==0 & dc == 0) continue; end
    %fprintf('dr = %d, dc = %d\n',dr,dc);

    % row and col of neighbors (rm,cm are the centers)
    rNbr = rm + dr; % Row of neighbor
    cNbr = cm + dc; % Col of neighbor
    
    % Make sure neighbors fall inside the slice
    indok = find(rNbr > 0 & rNbr <= nmrows & ...
		 cNbr > 0 & cNbr <= nmcols);
    rNbr = rNbr(indok);
    cNbr = cNbr(indok);

    % Keep track of the centers
    rm2 = rm(indok);
    cm2 = cm(indok);
    
    % Index into the slice of neighbors 
    indNbr = sub2ind(szm,rNbr,cNbr); 
    
    % Make sure neighbors fall inside the mask
    indok = find(mask(indNbr));
    
    % Final row, cols of neighbors and corresponding centers
    rNbr = rNbr(indok);
    cNbr = cNbr(indok);
    rm2 = rm2(indok);
    cm2 = cm2(indok);

    % Slice indices of neighbors and corresponding centers
    indNbr    = sub2ind(szm,rNbr,cNbr); 
    indCenter = sub2ind(szm,rm2,cm2);

    % Now get row and col of T. 
    rT = maskindex(indCenter);
    cT = maskindex(indNbr);
    
    % Now get indices of T
    indT = sub2ind(szT,rT,cT);
    
    % Set to 1
    T(indT) = 1;
  end
end

% Make each row sum to 1
Tsum = sum(T,2);
T = T ./ repmat(Tsum,[1 nm]);

% Finally ...
T = eye(nm) - T;

return;

