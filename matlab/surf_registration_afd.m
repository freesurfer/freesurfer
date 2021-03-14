function [prv ] = surf_registration_adf(subject, th_pval)

%% Tests the overall surface based registration %%
%


%
% surf_registration_afd.m
%
% Original Author: Laurence Wastiaux
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



hemi = sprintf('lh');
nv = 163842;

% a- Ressample the thickness curv file to sphere %
sn=strread(subject, '%s', 'delimiter', '/');
subject_name=char(sn(length(sn)));

trglh_file=sprintf('~/lh.thickness.average7.mgh');
trgrh_file=sprintf('~/rh.thickness.average7.mgh');

if(exist(trglh_file))
    cmdrm=sprintf('rm %s', trglh_file);
    u=unix(cmdrm);
end 
if(exist(trgrh_file))
    cmdrm=sprintf('rm %s', trgrh_file);
    u=unix(cmdrm);
end
%cmd_lh=sprintf('~/work/dev/mri_surf2surf/mri_surf2surf --srcsubject %s --srcsurfval thickness --src_type curv --trgsubject average7 --trgsurfval %s --trg_type mgh --hemi lh --noreshape',subject_name, trglh_file);
%cmd_rh=sprintf('mri_surf2surf --srcsubject %s --srcsurfval thickness --src_type curv --trgsubject average7 --trgsurfval %s --trg_type mgh --hemi rh --noreshape',subject_name, trgrh_file);
cmd_lh=sprintf('~/work/dev/mris_average_curvature/mris_average_curvature -a 100 -o average7 thickness lh sphere.reg %s %s', subject_name, trglh_file);
%cmd_rh=sprintf('mris_average_curvature -a 100 -o average7 thickness rh sphere.reg %s %s', subject_name, trgrh_file);
u1=unix(cmd_lh);
%u2=unix(cmd_rh);

% b- Load the ressampled surface, should be a (nv x 1) vector %
if (exist(trglh_file)==0) 
    fprintf('could not create the ressampled surface\n');
else
    %[d]=load_mgh(trglh_file);
    d=read_curv(trglh_file);
    if(size(d,1)~= nv)
        msgerr=sprintf('Wrong ressampled surface');
        error(msgerr)
    end
end
disp (size(d'))
d0=d';

%% Load all the data from the training set : the first evs can be loaded 
%% directly (lh.EV90.mgh has been computed for to remove 90% of the variance)

%load('/space/okapi/3/data/laurence/ADF/surf_registration/thickness_avg100new.mat'); %loads Dsreg
%Dsubj=load_mgh('/space/okapi/3/data/laurence/ADF/surf_registration/lh.ThicknessMatrix.mgh');

% Note: if the input subject belongs to the Buckner data set, the residual
% variance will be smaller than expected since it is included in Dsubj...

% demean
%thdemean = Dsubj - repmat(mean(Dsubj,2),[1 nv]);

% Covariance 
%M = thdemean*thdemean';

% SVD
%[u s v]=fast_svd(thdemean, M);
%ds=diag(s);
%ds2=diag(s.^2);
%fcvs=cumsum(ds/sum(ds));
%fcvs=cumsum(ds2/sum(ds2));
%nkeep=max(find(fcvs<.90))
%v=v(:, 1:nkeep);
v=load_mgh('/space/okapi/3/data/laurence/ADF/surf_registration/lh.EV90.mgh');
d1=d0-(d0*v)*v';

prv= var(d1)/var(d0);
[pinf psup]=compute_pval(prv);
if (psup<th_pval)
     mess=sprintf('The residual variance is too high (%g) %g %g', prv, pinf,psup);
     disp(mess)
else
     mess2=sprintf('Residual variance OK (%g) %g %g', prv,  pinf,psup);
     disp(mess2)
end

% subfunction compute_pval() %
function [p_inf, p_sup]=compute_pval(val)
%stat_file='/space/okapi/3/data/laurence/ADF/surf_registration/lh.SurfRegistrationResidualVar.adf'; 
%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end
stat_file=strcat(fsafdDir, '/lh.SurfRegistrationResidualVar.adf');
fid=fopen(stat_file);
if(fid==-1)
    mess=sprintf('Could not find %s', stat_file);
    error(mess)
end
while(strfind(fgetl(fid), '#'))
    pos=ftell(fid);
end
fseek(fid, pos, 'bof');
r=(fscanf(fid, '%g'))'; % r or r': does not make any difference
r=r/100;
pas=0.05;
x=0:pas:1;
[h] = hist(r,x);
p = h/sum(h);
dinf=find(x<=val);
dsup=find(x>val);
xinf=x(dinf);
xsup=x(dsup);
pinf=p(1:length(xinf));
psup=p(length(x)-length(xsup)+1:end);
if (val>=0 & length(xinf) >1 & length(xsup) >1)
    p_inf=trapz(xinf,pinf)/pas;
    p_sup=trapz(xsup,psup)/pas;
elseif (val>=0 & (length(xinf)<2 |  length(xsup)<2))
    pas2=pas/10;
    x2=0:pas2:1;
    [h2] = hist(r,x2);
    p2 = h2/sum(h2);
    dinf2=find(x2<=val);
    dsup2=find(x2>val);
    xinf2=x2(dinf2);
    xsup2=x2(dsup2);
    pinf2=p2(1:length(xinf2));
    psup2=p2(length(x2)-length(xsup2)+1:end);
    if(length(xinf2)>1)
        p_inf=trapz(xinf2,pinf2)/pas2;
    else
        p_inf=0;
    end
    if(length(xsup2)>1)
        p_sup=trapz(xsup2,psup2)/pas2;
    else
        p_sup=0;
    end  
else
    p_inf=0;
    p_sup=0;
end


