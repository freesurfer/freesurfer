function z = krfunc_lsqr(w,rimgsize,kimgsize,coilprof,transpstring)
% z = krfunc_lsqr(w,rimgsize,kimgsize,coilprof,transpstring)

z = [];
if(nargin < 4 | nargin > 5)
  fprintf('z = krfunc_lsqr(w,rimgsize,kimgsize,coilprof,transpstring)\n');
  error('');
  return;
end

if(nargin == 4)
  % No transpose, then run the forward model
  % where w = reconned image, and z = kspace image
  %fprintf('kimg:\n');
  z = kimgfunc(w,rimgsize,kimgsize,coilprof);
else
  % Transpose, then run the backprojection model
  % where w = kspace image, and z = reconned image
  %fprintf('rimg:\n');
  z = rimgfunc(w,rimgsize,kimgsize,coilprof);
end

return;
  
  





