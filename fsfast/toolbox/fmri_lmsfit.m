function h = fmri_lmsfit(Ch,ss)
%
% h = fmri_lmsfit(Ch,ss)
%
% Ch: nTotEst x nTotEst
% ss: nRows x nCols x nTotEst
% h:  nRows x nCols x nTotEst
%
% $Id: fmri_lmsfit.m,v 1.1 2003/03/04 20:47:40 greve Exp $

if(nargin ~= 2)
  msg = 'USAGE: fmri_lmsfit(Ch,ss)';
  qoe(msg);
  error(msg);
end

nTotEst = size(Ch,1);

szs = size(ss);
ns = length(szs);

% Determine nRows, nCols, and nTotEst %
if(ns == 2)
  nRows = size(ss,1);
  nCols = 1;
  nTotEst = size(ss,2);
else
  nRows = size(ss,1);
  nCols = size(ss,2);
  nTotEst = size(ss,3);
end

% Compute number of voxels %
nV = nRows*nCols; 

% Check dimensions of Ch and ss %%
if(size(Ch,1) ~= nTotEst | size(Ch,2) ~= nTotEst)
  msg = 'Ch and ss dimensions are inconsistent';
  qoe(msg);
  error(msg);
end

% Reshape ss for matrix mult %
ss = reshape(ss, [nV nTotEst])'; %'

% Heres where the action is %
h = Ch * ss;

% Reshape the output to the proper dim %
h = reshape(h', [nRows nCols nTotEst]); %'


return;

