function npernull = fast_npernull(nevents,nnullclocks)
%
% npernull = fast_npernull(nevents,nnullclocks)
%
% Given the total number of null clocks required (nnullclocks) 
% and the number of events, returns the number of null clocks
% that should occur before each event. One can think of it
% as inserting non-null events between each null event, so there
% needs to be one more null events than non-null events.
%
% For example, npernull(m) indicates the number of null clocks
% to place BEFORE the mth event. The last npernull is the number 
% of null clocks to place AFTER the last event.
%
% Sum(npernull) = nnullclocks. Length(npernull) = nevents + 1.
%
% Does not attempt to optimize the distribution of the lengths
% or to optimize the sequence of the lengths so that they are
% well counter-balanced. The distribution ends up being roughly
% exponential.

npernull = [];

if(nargin ~= 2)
  msg = 'USAGE: npernull = fast_npernull(nevents,nnullclocks)';
  fprintf('%s\n',msg);
  return;
end

% Randomly assign each null clock to a stimulus
r = ceil((nevents+1)*rand(nnullclocks,1));

% Count how many nulls were grouped with each stim
for n = 1:nevents+1,
  npernull(n) = length(find(r==n));
end

return;
