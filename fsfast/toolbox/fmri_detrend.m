function [ypost, coeff, hest] = fmri_detrend(ypre,X,Order,TPExclude)
%
% [ypost, coeff, hest] = fmri_detrend(ypre,X,Order,TPExclude)
%
% Removes all trends of order 0 (mean) to order Order-1 while
% simultanesouly fitting for the hemodynamic response. Detrending
% is done on a run-by-run basis.  If you do not want the HDR to
% be fit, pass X=[].  hest will be [] in that case.
%
% [ypost, coeff, hest] = fmri_detrend(ypre,X,Order)
%
% ypre:  raw fMRI slices (nRows x nCols x nTP x nRuns)
% X:     Stim Conv Mtx (nTP x nTotEst x nRuns)
% Order: number of trend components to remove (scalar-int)
% 
% ypost:   detrended fMRI slices (nRows x nCols x nTP x nRuns)
% coeff:   trend coefficients (nRows x nCols x Order x nRuns)
% hest:    HDR estimates for each run (nRows x nCols x Nch x nRuns)
%
%


%
% fmri_detrend.m
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

if(nargin ~= 4)
  msg = 'USAGE: [ypost, coeff] = fmri_detrend(ypre,X,Order,TPExclude)';
  qoe(msg);  error(msg);
end

if(Order == 0)
  fprintf(1,'INFO: Detrend Order = 0, no action taken\n');
  ypost = ypre;
  coeff = 0;
  hest = 0;
  return;
end

szy = size(ypre);
ny = length(szy);

nRows = size(ypre,1);
nCols = size(ypre,2);
nTP   = size(ypre,3);
nRuns = size(ypre,4);

if(~isempty(X))
  if(size(X,1) ~= nTP | size(X,3) ~= nRuns)
    msg = 'ypre and X dimensions are inconsistent';
    qoe(msg);
    error(msg);
  end
end

hest = [];

nTotEst = size(X,2);

% Compute number of voxels %
nV = nRows*nCols; 

ypre = reshape(ypre, [nV nTP nRuns]);
ypre = permute(ypre, [2 1 3]);

% Construct columns for detrending %
t = [0:nTP-1]'; %'
f = [];
for n = 0:Order-1,
  f = [f t.^n];
end

ypost = zeros(size(ypre));

for r = 1:nRuns,
  % fprintf(1,'   --- Detrending Run %d ---\n',r);

  y = ypre(:,:,r);          % get the rth functional slice

  if(~isempty(X)) SCM = [X(:,:,r) f]; % construct an SCM with detrend columns
  else            SCM = f;            % use only detrend columns
  end

  % Exclude specified data points %%
  iTPExclude = find(TPExclude(:,r)==1);
  if(~isempty(iTPExclude))
     nExl = length(iTPExclude);
     %fprintf(1,'   Run %d: Excluding %d Data Points:\n',r,nExl);
     %fprintf(1,'   ');
     %fprintf(1,'%d ',iTPExclude);
     %fprintf(1,'\n',iTPExclude);
     y(iTPExclude,:) = zeros(nExl,nV);
     SCM(iTPExclude,:)    = zeros(nExl,size(SCM,2));
  end

  %fprintf(1,'   Computing h\n');
  h = inv(SCM'*SCM)*SCM'*y; % estimate trend and HDR simultaneously 

  %fprintf(1,'   Extracting coeff\n');
  c = h([nTotEst+1:nTotEst+Order],:); % extract trend coeff

  %fprintf(1,'   Removing Trend\n');
  ypost(:,:,r) = y - f*c;   % remove trend
  %fprintf('     Done removing trend\n');
   if(~isempty(iTPExclude))
     ypost(iTPExclude,:,r) = zeros(nExl,nV);
   end

  %fprintf(1,'   Copying coeffs\n');
  if(nargout > 1)  coeff(:,:,r) = c; end
  %fprintf(1,'   Copying hest\n');
  if(nargout > 2 & ~isempty(X))  hest(:,:,r)  = h([1:nTotEst],:); end
end

clear ypre;
ypost = permute(ypost, [2 1 3]);
ypost = reshape(ypost, [nRows nCols nTP nRuns]);

if(nargout > 1)  
  coeff = permute(coeff, [2 1 3]);
  coeff = reshape(coeff, [nRows nCols Order nRuns]);
end

if(nargout > 2)  
  hest = permute(hest, [2 1 3]);
  hest = reshape(hest, [nRows nCols nTotEst nRuns]);
end

return;
