function XtX = fmri_idealXtX(nTypesPerSeq, Nh, FixationId)
%
% XtX = fmri_idealXtX(nTypesPerSeq, Nh, FixationId)
%
% Computes the ideal (X'*X) where X is the stimulus     %'
% convolution matrix.
% Note that the ideal normalized matrix is independent of the actual
% number of stimuli in the sequence.  Rather, it is dependent upon
% the probability of the stimuli.
%
% $Id: fmri_idealxtx.m,v 1.1 2003/03/04 20:47:39 greve Exp $


% If there is a fixation stim, remove fixation parameters
% and store in variable stps (Stim Types Per Run).
if(FixationId > -1)
  nStimTypes = length(nTypesPerSeq) - 1;
  indx = find([0:nStimTypes] ~= FixationId);
  stps = nTypesPerSeq(indx);
else % no Fixation
  nStimTypes = length(nTypesPerSeq);
  stps = nTypesPerSeq;
end

nStimTot = sum(stps);    % total number of Stim per Sequence
PrStim = stps/nStimTot;  % probability of each Stim Type

nXtX = nStimTypes * Nh;

XtX = zeros(nXtX,nXtX);

for rB=1:nStimTypes,
  for cB=1:nStimTypes,
    for rE=1:Nh,
      for cE=1:Nh,

        r = (rB-1)*Nh + rE;
        c = (cB-1)*Nh + cE;

        if(rB==cB) % diagonal block
          if(rE==cE) % diagonal element
	    XtX(r,c) = PrStim(rB);
          else % off-diagonal
	    XtX(r,c) = PrStim(rB)*PrStim(rB);
          end

	else % off-diagonal block
          if(rE==cE) % diagonal element
	    XtX(r,c) = 0;
          else % off-diagonal
	    XtX(r,c) = PrStim(rB)*PrStim(cB);
          end
        end

      end % for cE %
    end % for rE %
  end % for cB %
end % for rB %

return;



