function [vSig, pSig] = fmri_randeff(Test, h_is, hs, RM, Ms)
%
% [vSig, pSig] = fmri_randeff(Test, h_is, hs, RM, Ms)
%
% Random effects model for inter-session statistical analysis.
%  Test - test type (t or F) - only t works right now
%  h_is - average of the HDRs across subjects at each voxel
%  hs   - HDRs for each subject.
%  RM   - restriction matrix
%  Ms   - inter-session averaging matrix (see fmri_avgmtx_is)
%
%
%


%
% fmri_randeff.m
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

if(nargin ~= 5)
 msg = 'USAGE: [vSig, pSig] = fmri_randeff(Test, h_is, hs, RM, Ms)';
 qoe(msg);error(msg);
end

%%%% -------  Determine the Test Type ------------ %%%%%%
if(strncmp('T',upper(Test),1))     TestId = 0; % t-Test
elseif(strncmp('F',upper(Test),1)) TestId = 1; % F-Test
else  error(sprintf('Test Type %s is unknown',Test));
end

%%% -- get the dimensions --- %%%
[nRows nCols Nch Ns] = size(hs);
nV = nRows * nCols;
DOF = Ns - 1;
if(DOF < 1)
  msg = 'Must have at lest two sessions';
  qoe(msg);error(msg);
end

J =  size(RM,1); % Rows of RM %

%%%% ----- reshape for easier processing ----- %%%%%
hs = reshape(hs, [nV Nch Ns]);
hs = permute(hs, [2 1 3]);
h_is = reshape(h_is, [nV Nch])'; %'

% Compute the difference from the average %
ds = hs - repmat(h_is, [1 1 Ns]);

q = sqrt(Ns)*RM * h_is;
if(TestId == 0) vSig = zeros(J,nV); % t-Test
else            vSig = zeros(1,nV); % F-Test
end

% Process each voxel %
for v = 1:nV,

  % Compute the average covariance matrix across session %
  Ch = 0;
  for s = 1:Ns
    Ch = Ch + ds(:,v,s) * ds(:,v,s)'; %'
  end
  Ch = Ch/DOF;

  %%% ---- Compute inv of DOF/Desgin covariance mtx --- %%%%
  RChRt   = RM * Ch * RM'; %'
  ind = find(RChRt==0);
  RChRt(ind) = .0000000000001;

  %%% --- Perform Tests --- %%%
  if(TestId == 0) % t-Test
    dRChRt = diag(RChRt);
    vSig(:,v) = q(:,v) ./ sqrt(dRChRt);    
  else  % F-Test
    vSig(v) = q(:,v)' * inv(RChRt) * q(:,v); %'
  end
end

if(TestId == 1) vSig = vSig/J; end% F-Test

if(nargout == 2)
  if(TestId == 0) 
    pSig = tTest(DOF,reshape1d(vSig),75); % dof>75 --> normal approx %
  else
    pSig = FTest(J,DOF,vSig);
  end

  pSig = reshape(pSig',[nRows nCols size(vSig,1)]);
end

%% Reshape into image dimensions %%
vSig = reshape(vSig',[nRows nCols size(vSig,1)]);

return;
