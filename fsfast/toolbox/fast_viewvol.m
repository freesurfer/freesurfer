function h = fast_viewvol(vol,volres,frameres)
% h = fast_viewvol(vol,volres,frameres)

h = [];

if(nargin < 1 | nargin > 3)
  fprintf('USAGE: h = fast_viewvol(vol,volres,frameres)\n');
  return;
end

if(nargin < 2)      volres = [1 1 1]; end
if(isempty(volres)) volres = [1 1 1]; end

if(nargin < 3)        frameres = 1; end
if(isempty(frameres)) frameres = 1; end

