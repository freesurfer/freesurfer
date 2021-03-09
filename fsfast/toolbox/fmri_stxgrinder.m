function [vSig, pSig, ces] = fmri_stxgrinder(Test,hAvg, VoxVar, Ch, DOF, RM, q)
%
% [vSig, <pSig>, <ces> ] = fmri_stxgrinder(Test, hAvg, VoxVar, Ch, DOF, RM, q)
%
%
%
%


%
% fmri_stxgrinder.m
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


if(nargin < 6)
 msg = ...
'USAGE: [vSig, pSig, ces] = fmri_stxgrinder(Test, hAvg, VoxVar, Ch,DOF, RM , <q>)';
 qoe(msg);error(msg);
end

if(size(RM,2) ~= size(hAvg,3))
  fprintf('ERROR: dimension mismatch between the contrast matrix\n');
  fprintf('and the input averages. This may have happened if \n');
  fprintf('you redefined an analysis without redefining the \n');
  fprintf('contrast. Try re-running mkcontrast.\n');
  return;
end

%%%% -------  Determine the Test Type ------------ %%%%%%
if(strncmp('T',upper(Test),1))     TestId = 0; % t-Test
elseif(strncmp('F',upper(Test),1)) TestId = 1; % F-Test
else  error(sprintf('Test Type %s is unknown',Test));
end

%%% -- get the dimensions --- %%%
nRows = size(hAvg,1);
nCols = size(hAvg,2);
nT    = size(hAvg,3);
nV = nRows * nCols;

%% ---- Subtract the Ideal, if needed ----- %%%%%
if( nargin == 6 )
  h = hAvg;
else
  if(size(q,1) ~= size(hAvg,3))
    error('hAvg and q must have the same number of rows');
  end
  tmp = reshape(repmat(q,[nV 1]), [nRows nCols nT]);
  h = hAvg - tmp;
  clear tmp;
end

ind0 = find(VoxVar==0);
l0 = length(ind0);
fprintf(1,'  nVoxels with VoxVar=0: %3d\n',l0);
if(l0 == prod(size(VoxVar)))
  fprintf(1,'INFO: All voxels are zero\n');
  if(TestId == 0)
    sz = [nRows nCols size(RM,1)];
  else
    sz = [nRows nCols 1];
  end
  pSig = ones(sz);
  vSig = zeros(sz);
  ces =  zeros(sz);
  return;
end
if(l0 ~= 0)
  indnz = find(VoxVar ~= 0);
  VoxVar_min = min(reshape1d(VoxVar(indnz)));
  VoxVar(ind0) = VoxVar_min;
  fprintf(1,'  Setting zeros to min=%g\n',VoxVar_min);
end


%%%% ----- reshape for easier processing ----- %%%%%
h     = reshape(h,[nV nT])'; %'
VoxVarB = reshape(VoxVar,[nV 1])'; %'

%%% ---- Compute inv of DOF/Desgin covariance mtx --- %%%%
RChRt   = RM * Ch * RM'; %'
ind = find(RChRt==0);
RChRt(ind) = .0000000000001;

% Compute contrast effect size %
ces = RM*h;
nces = size(ces,1);
ces = reshape(ces',[nRows nCols nces]);

%%% --- Perform Tests --- %%%
if(TestId == 0) % t-Test
  fprintf('INFO: performing t-test\n');
  dRChRt = diag(RChRt);
  vSig = (RM * h) ./ sqrt(dRChRt * VoxVarB);

  pSig = tTest(DOF,reshape1d(vSig),300); % dof>300 --> normal approx %
  pSig = reshape(pSig,size(vSig));

else  % F-Test
  fprintf('INFO: performing F-test\n');
  if(strcmp('FM',upper(Test))) 
    dRChRt = diag(RChRt);
    vSig = (RM * h).^2 ./ (dRChRt * VoxVarB);
    J = 1;
  else
    RvvR = RM' * inv(RChRt) * RM; %'
    [J Nh] =  size(RM); % Rows of RM %
    if(Nh==1) vSig =  (((h' * RvvR)' .*h) ./ VoxVarB)/J;
    else      vSig =  ((sum((h' * RvvR)' .* h)) ./ VoxVarB)/J;
    end
  end

  pSig = FTest(J,DOF,reshape1d(vSig),1000); % 1500 = maxdof2
  pSig = reshape(pSig,size(vSig));
end

%% Reshape into image dimensions %%
vSig = reshape(vSig',[nRows nCols size(vSig,1)]); %'
pSig = reshape(pSig',[nRows nCols size(pSig,1)]); %'

return;
