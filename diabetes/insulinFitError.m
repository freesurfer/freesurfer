function [E, BG] = insulinFitError(k1, k3, k45, dG0, SC, SI, b1, b2, bTime, cTime, param)

sK = bsxfun(@plus, k3, 2*k45);
sD = sqrt(bsxfun(@minus, sK.^2, 4*k45.^2));
w1 = .5*(sK - sD);
a = .5*(sK + sD) ./ w1;
clear sK sD

cTimeD = bsxfun(@minus, param.time, cTime);
bolusTimeD = param.time - param.bolusTime;
bTimeD  = bsxfun(@minus, param.time, bTime);
pbTimeD = bsxfun(@minus, param.time, param.basalTime);

BG = bsxfun(@plus, param.BG(param.time==0) + dG0, bsxfun(@minus, param.C*bsxfun(@times, SC, bsxfun(@times, (1-exp(-bsxfun(@times, k1, cTimeD))), cTimeD>=0)), bsxfun(@times, SI, bsxfun(@plus, bsxfun(@plus, param.bolus*bsxfun(@times, 1-bsxfun(@rdivide, bsxfun(@times, a, exp(-bsxfun(@times, w1, bolusTimeD)))-exp(-bsxfun(@times, a.*w1, bolusTimeD)), a-1), bolusTimeD>=0), bsxfun(@times, param.b1 - b1, param.time)), bsxfun(@plus, bsxfun(@times, param.b2 - (param.b1 - b1), bsxfun(@times, bsxfun(@minus, pbTimeD, (1./w1).*(1+1./a)) + bsxfun(@rdivide, bsxfun(@times, a.^2, exp(-bsxfun(@times, w1, pbTimeD)))-exp(-bsxfun(@times, a.*w1, pbTimeD)), w1.*a.*(a-1)), pbTimeD>=0)), bsxfun(@times, bsxfun(@minus, -b2, param.b1 - b1), bsxfun(@times, bsxfun(@minus, bTimeD, (1./w1).*(1+1./a)) + bsxfun(@rdivide, bsxfun(@times, a.^2, exp(-bsxfun(@times, w1, bTimeD)))-exp(-bsxfun(@times, a.*w1, bTimeD)), w1.*a.*(a-1)), bTimeD>=0)))))));
E = mean(abs(bsxfun(@minus, BG, param.BG)), find(size(param.BG)>1));  % L1 error
%E = sqrt(mean(bsxfun(@minus, BG, param.BG).^2, find(size(param.BG)>1)));  % L2 error
% bsxfun(@plus, bsxfun(@plus, param.bolus*bsxfun(@times, 1-bsxfun(@rdivide, bsxfun(@times, a, exp(-bsxfun(@times, w1, bolusTimeD)))-exp(-bsxfun(@times, a.*w1, bolusTimeD)), a-1), bolusTimeD>=0), bsxfun(@times, param.b1 - b1, param.time)), bsxfun(@times, bsxfun(@minus, param.b2 - b2, param.b1 - b1), bsxfun(@times, bsxfun(@minus, bTimeD, (1./w1).*(1+1./a)) + bsxfun(@rdivide, bsxfun(@times, a.^2, exp(-bsxfun(@times, w1, bTimeD)))-exp(-bsxfun(@times, a.*w1, bTimeD)), w1.*a.*(a-1)), bTimeD>=0))) % metabolized insulin
%plot(param.time(:), BG(:), param.time(:), param.BG(:))
