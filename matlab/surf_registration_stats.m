function [Dsubj, nf, r, thdemean, M] = surf_registration_stats(dirname, arg2, arg3)


%
% surf_registration_stats.m
%
% Original Author: Laurence Wastiaux
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%



LoadTHMat=1;
D=[];
th=90;
if (nargin==3)
    D=arg2;
    th=arg3;
    LoadTHMat=0;
end
if (nargin==2)
    D=arg2;
    th=0.90;
    LoadTHMat=0;
end
files=dir(dirname);
nf=files;
length(files)
nv = 163842;
%D=zeros(length(files), nv);

%% 1- Load the thickness matrix (N_vertices x N_subjects) %%
if(LoadTHMat==1);
    for i_sub=5:300
        subject_name=files(i_sub).name;
        disp(subject_name)
        SubjectDir=strcat(dirname, '/', subject_name);
        SurfDir=strcat(SubjectDir, '/surf');
        
        
        % a- Ressample the thickness curv file to sphere %
        trglh_file=sprintf('~/lh.thickness.average7.mgh');
        trgrh_file=sprintf('~/rh.thickness.average7.mgh');
        
        if(exist(trglh_file))
            cmdrm=sprintf('rm %s', trglh_file);
            u=unix(cmdrm);
        end 
        if(exist(trglh_file))
            cmdrm=sprintf('rm %s', trglh_file);
            u=unix(cmdrm);
        end
        
        lhthickness_file=strcat(SurfDir, '/lh.thickness');
        rhthickness_file=strcat(SurfDir, '/rh.thickness');
        if (isempty(SurfDir) | ~exist(lhthickness_file) | ~exist(rhthickness_file) )
            fprintf('thickness files do not exist\n');
            i_sub=i_sub+1;
        else
            %cmd_lh=sprintf('mri_surf2surf --srcsubject %s --srcsurfval thickness --src_type curv --trgsubject average7 --trgsurfval %s --trg_type mgh --hemi lh --nsmooth-in 100 --noreshape',subject_name, trglh_file);
            %cmd_rh=sprintf('mri_surf2surf --srcsubject %s --srcsurfval thickness --src_type curv --trgsubject average7 --trgsurfval %s --trg_type mgh --hemi rh --nsmooth-in 100 --noreshape',subject_name, trgrh_file);
            cmd_lh=sprintf('~/work/dev/mris_average_curvature/mris_average_curvature -a 100 -o average7 thickness lh sphere.reg %s %s', subject_name, trglh_file);
            u1=unix(cmd_lh);
            %u2=unix(cmd_rh);
            
            % b- Load the ressampled surface %
            if (exist(trglh_file)==0)
                fprintf('could not create the ressampled surface\n');
                i_sub=i_sub+1;
            else
                %[vol]=load_mgh(trglh_file);
                vol=read_curv(trglh_file);
                sz=size(vol)
                if(sz(1)~= nv)
                    msgerr=sprintf('Wrong ressampled surface');
                    error(msgerr)
                end
                D= [D ; vol'];
                nf(size(D,1)).name=subject_name;
            end
        end
    end
end
Dsubj=D;
%save('/space/okapi/3/data/laurence/surf_registration/thickness_avg100.mat', 'Dsubj');
%save('/space/okapi/3/data/laurence/surf_registration/Subj_Names.mat', 'nf');
sd=size(D);

nsubjs=sd(1);
nns=1:nsubjs;

%% 2- compute stats %%

% demean
thdemean = Dsubj - repmat(mean(Dsubj,2),[1 nv]);

% Covariance 
M = thdemean*thdemean';

% % Jackknife
for jksubj = 1:nsubjs
  fprintf('%d *\n',jksubj);
  indjk = find(nns ~= jksubj);

  % SVD 
  [u s v] = fast_svd(thdemean(indjk,:),M(indjk,indjk)); 
  fcvs(:,jksubj) = cumsum((diag(s.^2))/(sum(diag(s.^2))));
  
  % keep enough to explain 85% of the data. The fcvs does not
  % change that much regardless of which data set is jked
  nkeep(jksubj) = max(find(fcvs(:,jksubj)<th));
  v = v(:,1:nkeep(jksubj));
  
  % Extract the jkth subjects thickness and compute var
  thsubj = thdemean(jksubj,:);
  %thsubj = randn(size(thsubj));
  d0(jksubj) = var(thsubj);

  % Project out the first nkeep EVs and recompute var
  thsubjproj = thsubj - (thsubj*v)*v';
  d1(jksubj) = var(thsubjproj);
  
  fprintf('  nkeep = %3d, d0=%g, d1=%g, r=%g, rwhite = %g\n',...
	  nkeep(jksubj),d0(jksubj),d1(jksubj),d1(jksubj)/d0(jksubj),1-nkeep(jksubj)/size(thdemean,1));

end

% Compute percent residual variance
r = 100*d1./d0;
%save('/space/okapi/3/data/laurence/ADF/surf_registration/residual_var80.mat', 'r');
% Plot the histogram
% indr = find(r<25);
% rr = r(indr);
 [h x] = hist(r);
 p = h/sum(h);
 ph = plot(x,p,'r-');
% set(ph,'Linewidth',3); 
% xlabel('Percent Residual Variance');
% title('Distribution of Percent Residual Variance in Thickness')
% set(gcf,'color',[1 1 1]);

