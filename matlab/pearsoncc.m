function rho = pearsoncc(x,y)
% rho = pearsoncc(x,y)
% x and y should be column vectors
% If x and y have multiple columns, then
%   rho is computed for each one

rho = [];
if(size(x,1) ~= size(y,1))
  fprintf('ERROR: x and y have different number of rows\n');
  return;
end
if(size(x,2) ~= size(y,2))
  fprintf('ERROR: x and y have different number of columns\n');
  return;
end

% Remove the mean
xmn = mean(x);
dx = x - repmat(xmn,[size(x,1) 1]);
ymn = mean(y);
dy = y - repmat(ymn,[size(y,1) 1]);

rho = sum(dx.*dy)./sqrt(sum(dx.*dx).*sum(dy.*dy));

return;

