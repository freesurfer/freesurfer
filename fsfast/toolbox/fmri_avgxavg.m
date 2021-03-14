function [havgavg, hstdavg, DOF, Ms] =fmri_avgxavg(effect, havg, hstd, xtx)
%
% [havgavg, hstdavg, DOF, Ms] = fmri_avgxavg(effect, havg, hstd)
% [havgavg, hstdavg, DOF, Ms] = fmri_avgxavg(effect, havg, hstd, xtx)
%


%
% fmri_avgxavg.m
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

if(nargin ~= 3 & nargin ~= 4)
   msg = 'USAGE: fmri_avgxavg(effect, havg, hstd, <xtx>)';
   qoe(msg);error(msg);
end

if( isempty( strmatch(lower(effect),{'random','fixed'},'exact') ) )
   msg = sprintf('%s Effect not recognized (use random or fixed)',effect);
   qoe(msg);error(msg);
end

if(nargin == 3) weighting = 'subject';
else            weighting = 'event';
end 

if( size(havg) ~= size(hstd) )
  msg = 'havg and hstd have incompatible dimensions';
  qoe(msg);error(msg);
end

[nRows nCols Nch Ns] = size(havg);
Nv = nRows*nCols;

if(strcmp(weighting,'event'))
  if(size(xtx,1) ~= Nch | size(xtx,3) ~= Ns)
    msg = 'dimension of xtx is not compatible with havg and hstd';
    qoe(msg);error(msg);
  end
end

effect = lower(effect);
DOF = Ns-1;

%%%%%%%%%%% Determine Weighting Matrix %%%%%%%%%%%%%%%
if(strcmp(weighting,'subject')) % subject weighted %
   Ms = eye(Nch)/Ns;
   Ms = repmat(Ms, [1 1 Ns]);
else   %'event'%           % event-weighted %
   dxtxsum = 0;
   Ms = zeros(Nch,Nch,Ns);
   for s=1:Ns, 
     dxtx(:,s) = diag(xtx(:,:,s));
     dxtxsum = dxtxsum + dxtx(:,s);
   end
   for s=1:Ns, 
     v = diag(xtx(:,:,s)) ./ dxtxsum;
     Ms(:,:,s) = diag(v);
   end
end

%%%%%%%%%%% Average %%%%%%%%%%%%%%%%%%%%%%%%%%%%
havgavg = zeros(Nch,Nv);
hvaravg = zeros(Nch,Nv);
havg = permute(havg, [3 1 2 4]);
havg = reshape(havg, [Nch Nv Ns]);

%% Compute pooled averages, given the weighting matrix %%%%
for s=1:Ns,
   havgavg = havgavg + Ms(:,:,s) * havg(:,:,s);
end

%%%Compute pooled variances according to model %%%%%
if(strcmp(effect,'random'))
  havgerr = havg - repmat(havgavg, [1 1 Ns]);
  hstdavg = std(havgerr,[],3);

else % 'fixed' effects %
  hstd = permute(hstd, [3 1 2 4]);
  hstd = reshape(hstd, [Nch Nv Ns]);
  hvar = hstd.^2;
  clear hstd;

  for s=1:Ns,
     hvaravg = hvaravg + Ms(:,:,s) * hvar(:,:,s);
  end

  hstdavg = sqrt(hvaravg);
  clear hvaravg hvar;

end

havgavg = reshape(havgavg',[nRows nCols Nch]); %'
hstdavg = reshape(hstdavg',[nRows nCols Nch]); %'




return
