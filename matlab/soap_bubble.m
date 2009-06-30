function Ismooth = soap_bubble(I, ctrl, niter)
% Ismooth = soap_bubble(I, ctrl, niter)
% or
% Ismooth = soap_bubble(I, niter)
% in which case the nonzero entries of I will be taken to be
% the control points

if (nargin < 3)
  niter = ctrl ;
  ind = find(I ~= 0) ;
  ctrl = zeros(size(I)) ;
  ctrl(ind) = ones(size(ind)) ;
end
  
Ismooth = I ;

ind = find(ctrl ~= 0) ;
filter = fspecial('average') ;
  
for  n=1:niter
  Ismooth = filter2(filter, Ismooth, 'same') ;
  Ismooth(ind) = I(ind) ;
end



