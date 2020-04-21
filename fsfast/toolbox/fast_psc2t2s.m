function T2sAct = fast_psc2t2s(PctSigChange,TE,T2sBaseline)
% T2sAct = fast_psc2t2s(PctSigChange,TE,T2sBaseline)
%
% Computes the T2-star during activation given the percent signal
% change of the activation (relatitve to baseline), the echo time,
% and the T2-star at baseline.
%

T2sAct = [];
if(nargin ~= 3)
  fprintf('T2sAct = fast_psc2t2s(PctSigChange,TE,T2sBaseline)\n');
  return;
end

RelSigChange = PctSigChange/100;
R2sBaseline = 1/T2sBaseline;
dR2s = -log(RelSigChange)/TE;
R2sAct = R2sBaseline + dR2s;
T2sAct = 1/R2sAct;

return;

