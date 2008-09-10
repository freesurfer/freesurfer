function indtrig = triggerfinder(trigger)
% indtrig = triggerfinder(trigger)
% 
% Finds the index of the rising edge of the trigger. Trigger should
% be a (mostly) binary waveform. The threshold is set to half way
% between the min and max.
%
% $Id: triggerfinder.m,v 1.1 2008/09/10 18:15:50 greve Exp $

indtrig = [];
if(nargin ~= 2)
  fprintf('indtrig = triggerfinder(trigger)\n');
  return
end

thresh = (max(trigger)-min(trigger))/2.0;

trigger = trigger > thresh;
dtrigger = diff(trigger);
indtrigger = find(dtrigger == 1) + 1;

return;







