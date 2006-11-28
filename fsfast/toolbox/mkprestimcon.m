function C = mkprestimcon(sxa,wcond)
% C = mkprestimcon(sxa,wcond)

C = [];
if(nargin ~= 2)
  fprintf('C = mkprestimcon(sxa,wcond)\n');
  return;
end

nPre = round(sxa.TPreStim/sxa.TER);
nFIR = sxa.Nh;
Nc = sxa.Nc-1;
wcond = wcond(:)';

if(length(wcond) ~= Nc)
  fprintf('ERROR: wcond has wrong number of items\n');
  return;
end

a = zeros(nFIR);
a(:,1:nPre) = -1/nPre;
b = a + eye(nFIR);
C = kron(wcond,b); 

return;
