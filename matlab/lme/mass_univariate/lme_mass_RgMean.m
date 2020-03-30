function [RgMeans,nRg] = lme_mass_RgMean(Rgs,Data)
% [RgMeans,nRg] = lme_mass_RgMean(Rgs,Data)
% 
% Mean of the data for each region.
%
% Input
% Rgs: 1 x nv segmentation vector containing a region number assigned 
% to each vertice along the surface.
% Data: Any data (eg. a matrix whose colums are estimators of the covariance
% components at each location.
%
% Output
% RgMeans: Matrix assigning the mean of Data within a region to each vertex
% belonging to that region.
% nRg: Number of regions in Rgs.
%
% Original Author: Jorge Luis Bernal Rusiel 
%
if nargin < 2
    error('Too few inputs');
end;
[npar,nv] = size(Data);
if length(Rgs) ~= nv
    error('Both inputs must have the same number of colums');
end;
RgMeans = zeros(npar,nv);
Rgnums = unique(Rgs);
nRg = length(Rgnums);
for i=1:nRg
    Rgvtxs = find(Rgs == Rgnums(i));
    mRgData = mean(Data(:,Rgvtxs),2);
    RgMeans(:,Rgvtxs) = kron(ones(1,length(Rgvtxs)),mRgData);
end;
