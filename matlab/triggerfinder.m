function indtrigger = triggerfinder(trigger)
% indtrigger = triggerfinder(trigger)
% 
% Finds the index of the rising edge of the trigger. Trigger should
% be a (mostly) binary waveform. The threshold is set to half way
% between the min and max.
%

indtrig = [];
if(nargin ~= 1)
  fprintf('indtrigger = triggerfinder(trigger)\n');
  return
end

thresh = (max(trigger)-min(trigger))/2.0;

trigger = trigger > thresh;
dtrigger = diff(trigger);
indtrigger = find(dtrigger == 1) + 1;

return;







