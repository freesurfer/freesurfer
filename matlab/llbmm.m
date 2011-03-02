function ll = llbmm(nn,pA,pI,lambda,M)
% ll = llbmm(nn,pA,pI,lambda,<M>)
%
% Computes the log-likelihood of data (nn) given a binomial mixture
% model with parameters pA, pI, and lambda.
%
% nn is the histogram over all voxels of the number of times
% that a given a given number of positives was detected. The
% number of trials M is assumed to be the length of nn, unless
% specfied explicitly with M. 
%
% pA - probability of declaring a voxel to be active given that it
% is truly active (the TPR). Eg, .8
%
% pI - probability of declaring a voxel to be active given that it
% is truly inactive (the FPR). Eg, .01 (the voxel-wise threshold).
%
% lambda - proportion of truly active voxels out of all voxels
%
% Note that pA != 1 - pI.
%
% Based on Genovese, et al, 1997. Estimating Test-Retest
% Reliability in Functional MR Imaging I: Statistical Methodology. 
% MRM 38:497-507. 
%
% See the bottom of this file for testing code.
%


%
% llbmm.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.4 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


ll = [];
if(nargin < 4 | nargin > 5)
  fprintf('ll = llbmm(nn,pA,pI,lambda,<M>)\n');
  return;
end

% Number of trials
if(~exist('M','var')) M = []; end
if(isempty(M)) M = length(nn)-1; end

nn = nn(:); % assure a column vector
if(M ~= length(nn)-1)
  if(M < length(nn)-1)
    fprintf('ERROR: length of nn (%d) > M (%d)\n',...
	    length(nn),M);
    return;
  end
  % Pad nn with zeros
  nn = [nn; zeros(M-length(nn),1)];
end

k = [0:M]'; % All possible outcomes

a = lambda .* (pA.^k) .* ((1-pA).^(M-k)) ;
b = (1-lambda) .* (pI.^k) .* ((1-pI).^(M-k)) ;

ll = sum(nn .* log(a + b));

return;

%-------------------------------------------------------------------%
% Below is code to test that this function does produce a maximum at
% the expected parameter set.

Na =   100;  % Number of active voxels
Ni = 10000;  % Number of inactive voxels
N = Na+Ni; % Total number of voxels
lambda = Na/N; % ratio of active to total

pA = .2; % Probability that a truly active vox is detected (TPR)
pI = .1; % Probability that a truly inactive vox is detected (FPR)
% Note: pA+pI != 1

M = 40; % Number of trials

% Number of trials each of the Na voxels was delcared active. This
% will be a list of Na numbers, each number between 0 and M.
A = randb(pA,M,Na);

% Number of trials each of the Ni voxels was delcared active. This
% will be a list of Ni numbers, each number between 0 and M.
I = randb(pI,M,Ni);

% This is the synthesized "data", ie, all active and inactive voxels
% mixed together, which the value at each voxel is the number of times
% it was declared active over the M trials.
D = [A;I];

% Build a histogram of D
x = [0:M]; % This is a list of all possible values
h0 = hist(D,x); % Count
h = h0/N; % Probability

% Construct ideal PDF of the mixture under these conditions
pdfA = binomialpdf(x,M,pA);
pdfI = binomialpdf(x,M,pI);
pdf = lambda*pdfA + (1-lambda)*pdfI;

% Plot the ideal vs the actual
plot(x,pdf,'+-',x,h,'*-');

% This is a test in which the log-likelihood of the synthesized data
% is computed for various values of pA to test whether the max occurs
% at the ideal value.
pAlist = [.01:.01:.99];
clear ll;
for nth = 1:length(pAlist);
  ll(nth) = llbmm(h0,pAlist(nth),pI,lambda);
end
[m k] = max(ll);
pAmax = pAlist(k);
fprintf('pA = %g, max = %g\n',pA,pAmax);

% This is a test in which the log-likelihood of the synthesized data
% is computed for various values of pI to test whether the max occurs
% at the ideal value.
pIlist = [.01:.01:.99];
clear ll;
for nth = 1:length(pIlist);
  ll(nth) = llbmm(h0,pA,pIlist(nth),lambda);
end
[m k] = max(ll);
pImax = pIlist(k);
fprintf('pI = %g, max = %g\n',pI,pImax);

% This is a test in which the log-likelihood of the synthesized data
% is computed for various values of lambda to test whether the max occurs
% at the ideal value.
lambdalist = [.01:.01:.99];
clear ll;
for nth = 1:length(lambdalist);
  ll(nth) = llbmm(h0,pA,pI,lambdalist(nth));
end
[m k] = max(ll);
lambdamax = lambdalist(k);
fprintf('lambda = %g, max = %g\n',lambda,lambdamax);


% Search. Init is important because there is a symmetry between the
% parameters, ie, you get the same cost if you swap pA and pI and
% use 1-lambda.
params0 = [pA pI lambda];
initparams = [.5 .1 .5];
[optparams cost] = fminsearch('bmmcost',initparams,[],h0);
[initparams;params0;optparams]

[fpr,tpr,auc,fpr] = bmmroc(pA,pI,lambda,M);
plot(fpr,tpr)

