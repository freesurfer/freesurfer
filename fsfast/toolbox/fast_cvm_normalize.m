function [ocvm, acoravg, acorstd] = fast_cvm_normalize(icvm,nmax)
% ocvm = fast_cvm_normalize(icvm,nmax)
%
% Covariance matrix normalization.
%
% icvm - input covariance matrix
% ocvm - normalized covariance matrix
%
% Forces all the diagonals elements of ocvm to be equal. The value of the 
% output diagonal is equal to the mean of the input diagonal divided by
% the mean of the main diagonal (this forces the elements of the main 
% diagonal to be one). If nmax is set, then diagonals beyond nmax are 
% set to 0. Setting the diagonals to be equal enforces the assumption
% that the noise is stationary across observations.  
%
% Each diagonal is set to a constant equal to the mean of that diagonal
% divided by the mean of the main diagonal. This forces the main diagonal
% to be 1.
%


%
% fast_cvm_normalize.m
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

if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: ocvm = fast_cvm_stationary(icvm,nmax)';
  qoe(msg);error(msg);
end

if(size(icvm,1) ~= size(icvm,2))
  msg = 'Input cvm is not square';
  qoe(msg);error(msg);
end

ncvm = size(icvm,1);

% If nmax is not specified, use all diagonals %
if(nargin == 1) nmax = ncvm; end

% Make sure nmax is not too large %
if(nmax > ncvm)
  msg = sprintf('nmax = %d > ncvm = %d',nmax,ncvm);
  qoe(msg);error(msg);
end

% Get the mean of the main diagonal %
k0mean = mean(diag(icvm));

% Go through the off-center diagonals and compute the mean %
ocvm = zeros(ncvm);
for k = 1:nmax-1
  kdiag = diag(icvm,k)/k0mean;        % kth diagonal
  mnkdiag = mean(kdiag);              % mean of kth diagonal
  kdiag2 = mnkdiag*ones(size(kdiag)); % mean replicated into vector
  tmpcvm = diag(kdiag2,k);            % matrix with only kth diag set
  ocvm = ocvm + tmpcvm;               % accumulate
end

% make symmetric %
ocvm = ocvm + ocvm'; %'

% Set main diagonal to 1%
ocvm = ocvm + eye(ncvm);

return;
