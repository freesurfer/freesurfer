function [cc, yvar] = fast_corrcoef(y,lag,DOF)
% [cc, yvar] = fast_corrcoef(y,<lag>,<DOF>)
%
% y is nframes by ncols
% lag defaults to 1 if not present or null
% DOF defaults to nframes if not present or null
%
% cc is 1 by ncols set of normalized correlation coefficients
% yvar is 1 by ncols
%

cc = [];
yvar = [];

if(nargin < 1 | nargin > 3)
  fprintf('[cc, yvar] = fast_corrcoef(y,<lag>,<DOF>)\n');
  return; 
end

[nf nv] = size(y);

if(nargin < 2)   lag = 1; end
if(isempty(lag)) lag = 1; end
if(nargin ~= 3)  DOF = nf; end
if(isempty(DOF)) DOF = nf; end

if(DOF-lag <= 0)
  fprintf('ERROR: DOF-lag <= 0\n');
  return;
end

nn1 = 1:nf-lag;
nn2 = lag+1:nf;

yvar = sum(y.^2,1)/DOF;
ycvar = sum(y(nn1,:).*y(nn2,:),1)/(DOF-lag);
cc = ycvar./yvar;

return;
