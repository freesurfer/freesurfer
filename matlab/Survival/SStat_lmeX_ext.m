function [sID_ext,X_ext,d_ext,t_ext,Y_ext] = SStat_lmeX_ext(Rfx,nrfx,Bhat,sID,X,Scols,d,t)
% [sID_ext,X_ext,d_ext,t_ext,Y_ext] = SStat_lmeX_ext(Rfx,nrfx,Bhat,sID,X,Scols,d,t)
%
% This function is intended for the joint modelling of longitudinal and 
% time-to-event data. It takes as input the output from lme_mass_rfx plus 
% the time-to-event data and generates the requiered data for subsequent 
% application of the Extended Cox model (SStat_CoxExt). 
%
% Input
% Rfx: Estimated subject-especific random effects.
% nrfx: Number of random effects.
% Bhat: Population-level regression coefficient estimates.
% sID: Subjects' IDs (for each row of X).
% X: Longitudinal ordered design Matrix (according to time for each subject).
% Scols: Colums of X that are going to be considered as covariates for the
% subsequent extended Cox model. 
% d: Logical vector indicating time-invariant censorship status (1 if the 
% subject got the failure event during follow-up or 0 otherwise).
% t: Vector whose entries are the survival and censored times (ordered 
% according to X and time-invariant within each subject).
%
% Output
% sID_ext: Extended subjects' IDs.
% X_ext: Extended design matrix. Time-dependent covariates resulting from
% the interaction of the time-independent covariates in X_ext and t_ext can 
% be easily added.
% d_ext: Extended censorship status vector.
% t_ext: Extended survival time vector.

if nargin < 8 
    error('Too few inputs');   
end;
m = size(X,1);
if (length(sID)~=m) || (length(d)~=m) || (length(t)~=m) 
    error(['All, the design matrix X, data y, censorship status vector d, '...
        'time vector t and subject ID vector must have the same number of rows.']);
end;
[usID,ix_sID] = unique(sID,'stable');
display('Generating the extended Cox model data');
X_Surv = X(ix_sID,Scols);
t_Surv = t(ix_sID);
d_Surv = d(ix_sID);
[sID_ext,X_ext,d_ext,t_ext] = SStat_X_ext(usID,[Rfx X_Surv],d_Surv,t_Surv);
nv = size(Bhat,2);
Y_ext = zeros(size(X_ext,1),nv);
for i=1:nrfx
  Y_ext = Y_ext + (kron(ones(size(X_ext,1),1),Bhat(i,:)) + X_ext(:,i:nrfx:nv*nrfx-nrfx+i)).*kron(ones(1,nv),(t_ext.^(i-1)));
end;
X_ext = X_ext(:,nv*nrfx+1:end);
display('Done');
