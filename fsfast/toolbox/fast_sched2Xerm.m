function X = fast_sched2Xerm(tPres,ntrs,TR,W,ERM)
% X = fast_sched2Xerm(tPres,ntrs,TR,W,ERM)
%
% Computes the regressor matrix for a given stimulus
% schedule and Event Response Model.
%
% ERM - is an event-response model structure. The
% structure has (at least) the following fields:
%   name - string with the name of the model
%   params - vector with parameter list
%
% Possible ERMs and their parameters:
%  fir         - tprestim, ter, timewindow
%  gamma       - delay, dispersion, boxcarwidth
%  gamma+deriv - delay, dispersion, boxcarwidth
%
% See also:
%   fast_sched2Xfir, fast_sched2Xgamma, fast_sched2Xgammaderiv
%   


%
% fast_sched2Xerm.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

X = [];

if(nargin ~= 5)
  fprintf('X = fast_sched2Xerm(tPres,ntrs,TR,W,ERM)\n');
  return;
end

if(~isfield(ERM,'name'))
  fprintf('ERROR: ERM does not have a name field\n');
  return;
end

if(~isfield(ERM,'params'))
  fprintf('ERROR: ERM does not have a params field\n');
  return;
end

nparams = length(ERM.params);

switch(lower(ERM.name))

   case 'fir',
     if(nparams < 3)     
       fprintf('ERROR: ERM fir requires 3 parameters (found %d)\n',nparams);
       return;
     end
     TPreStim   = ERM.params(1);
     TER        = ERM.params(2);
     TimeWindow = ERM.params(3);
     X = fast_sched2Xfir(tPres,ntrs,TR,TER,TPreStim,TimeWindow,W);

   case 'gamma',
     if(nparams < 3)     
       fprintf('ERROR: ERM gamma requires 3 parameters (found %d)\n',nparams);
       return;
     end
     Delay       = ERM.params(1);
     Dispersion  = ERM.params(2);
     BoxCarWidth = ERM.params(3);
     X = fast_sched2Xgamma(tPres,ntrs,TR,Delay,Dispersion,BoxCarWidth,W);

   case 'gamma+deriv',
     if(nparams < 3)     
       fprintf('ERROR: ERM gamma+deriv requires 3 parameters (found %d)\n',nparams);
       return;
     end
     Delay       = ERM.params(1);
     Dispersion  = ERM.params(2);
     BoxCarWidth = ERM.params(3);
     X = fast_sched2Xgammaderiv(tPres,ntrs,TR,Delay,Dispersion,BoxCarWidth,W);

   otherwise,
     fprintf('ERROR: Event Response Model %s unrecognized\n',ERM.name);
       return;

end


return;
