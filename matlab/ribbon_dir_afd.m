function [D, Isubj]=ribbon_dir_adf(dirname, th_pval)

%      For each subject in the directory "dirname":
% Computes the Dice coefficients measuring the overlap of the
% Cortical Ribbon volume computed 
%         1- from the subcortical labeling   
%         2- as the space between the white and the pial surface 
%


%
% ribbon_dir_afd.m
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




if (nargin<1 | nargin>2)
    msg=sprintf('USAGE: [D]=ribbon_dir_adf(dirname, th_pval)');
    disp(msg)
end

%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end

% Load stats from the Buckner data set %
%stat_file='/space/okapi/3/data/laurence/ADF/ribbon/RibbonPlacementDice.adf';
stat_file=strcat(fsafdDir, '/RibbonPlacementDice.adf');
fid=fopen(stat_file);
if(fid==-1)
    mess=sprintf('Could not find %s', stat_file);
    error(mess)
end
while(strfind(fgetl(fid), '#'))
    pos=ftell(fid);
end
fseek(fid, pos, 'bof');
DD=fscanf(fid, '%g');

files=dir(dirname);
flag=1;
D=[];

for i=5:length(files)
    flag=1;
    SubjectDir=strcat(dirname,'/',files(i).name);
    i_subj=i;
    disp(SubjectDir)
    create_ribbon=0;
    
    %%  1- Load the aseg volume    %%
    %% --------------------------  %%
    
    CorDir1=strcat(SubjectDir,'/mri/aseg/');
    d1=dir(CorDir1);

    if (length(d1)<3 | ( length(strfind(files(i).name,'0'))==0 &&(length(strfind(files(i).name,'1'))==0 )))
        mess=sprintf('Cannot find the coronals of the volume aseg for the subject %s',files(i).name);
        disp(mess)
        flag=0;
    else
        [vol1 M]=load_cor2(SubjectDir,'aseg');
    
    end

    %%  2- Look for ribbon mask or create it if not present %%
    %% ---------------------------------------------------  %%
    
    % a- Check if the left and right surfaces exist %%
    SurfDir=strcat(SubjectDir, '/surf');
    surf=dir(SurfDir);
    if(length(surf)<=2)
        flag=0;
        mess=sprintf('Directory %s is empty', SurfDir);
        disp(mess)
    end
    if(length(surf)>2)
        if (~(exist(strcat(SurfDir, '/rh.ribbon.mgh'))) | ~(exist(strcat(SurfDir, '/lh.ribbon.mgh'))) )
        mess=sprintf('Cannot find the cortical ribbon for subject %s... try to create it', files(i).name);
        disp(mess)
        create_ribbon=1;
        end
    end
    % a- Check if the left and right surfaces exist %%
    if(length(surf)>2)
        if ((exist(strcat(SurfDir, '/rh.white'))==0) | (exist(strcat(SurfDir, '/lh.white'))==0) | (exist(strcat(SurfDir, '/rh.thickness'))==0) | (exist(strcat(SurfDir, '/lh.thickness'))==0))
            flag=0;
            mess=sprintf('Cannot find the white surface for subject %s', files(i).name);
            disp(mess)
        end
    end
     
    if(create_ribbon==1 & flag~=0) %Create cortical ribbon
        flag2=1;
        % b- Generate a volume registration file 
        a=strread(SubjectDir, '%s', 'delimiter', '/');
        subject_name=char(a(length(a)));
        reg_file=sprintf('~/reg.dat');
        fid=fopen(reg_file, 'w');
        fprintf(fid, '%s\n', subject_name);
        fprintf(fid, '1\n1\n1\n1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1\n%s\n', 'round');
        fclose(fid);
        
        % c- Fill ribbon in the left & right hemispheres
        template_vol=strcat(SubjectDir, '/mri/orig/');
        lh_ribbon=sprintf('~/ribbon_lh.mgz');
        rh_ribbon=sprintf('~/ribbon_rh.mgz');
        
        if( exist(lh_ribbon)~=0)
            cmdrm1=sprintf('rm %s', lh_ribbon);
            u1=unix(cmdrm1);
        end
        if( exist(rh_ribbon)~=0)
            cmdrm2=sprintf('rm %s', rh_ribbon);
            u1=unix(cmdrm2);
        end
        
        if(length(a)>1)
            sd=a(1);
            for i=1:(length(a)-2)
                sd=strcat(sd,'/', a(i+1));
            end
            SubjectsDir=char(sd);
        else
            SubjectsDir = deblank(getenv('SUBJECTS_DIR'));
            if(isempty(SubjectsDir))
                msg = 'Cannot find SUBJECTS_DIR environment variable';
                disp(msg);
                flag2=0;
            end
        end
    
        if(flag2~=0)
            cmd1=sprintf('~/work/dev/mri_surf2vol/mri_surf2vol --mkmask --hemi lh --fillribbon --template %s --outvol %s --volreg %s --sd %s', template_vol, lh_ribbon, reg_file, SubjectsDir);
            cmd2=sprintf('~/work/dev/mri_surf2vol/mri_surf2vol --mkmask --hemi rh --fillribbon --template %s --outvol %s --volreg %s --sd %s', template_vol, rh_ribbon, reg_file, SubjectsDir);
            
            c1=unix(cmd1);
            c2=unix(cmd2);
        end
    else  %if creat_ribbon==0
        if(exist(strcat(SubjectDir,'/surf/lh.ribbon.mgh')) & exist(strcat(SubjectDir,'/surf/rh.ribbon.mgh')) )
            lh_ribbon=strcat(SubjectDir,'/surf/lh.ribbon.mgh');
            rh_ribbon=strcat(SubjectDir,'/surf/rh.ribbon.mgh');
        end
        if(exist(strcat(SubjectDir,'/surf/lh.ribbon.mgz')) & exist(strcat(SubjectDir,'/surf/rh.ribbon.mgz')) )
            lh_ribbon=strcat(SubjectDir,'/surf/lh.ribbon.mgz');
            rh_ribbon=strcat(SubjectDir,'/surf/rh.ribbon.mgz');
        end
    end

    if(flag~=0)
        % d- Load the right & left cortical ribbon
        if( (exist(lh_ribbon)==0) | (exist(rh_ribbon)==0) )
            mess=sprintf('Could not generate cortical ribbon volume');
            disp(mess)
            i=i+1;
        else
            [vol_lhr]=load_mgh(lh_ribbon);
            [vol_rhr]=load_mgh(rh_ribbon);
            if(length(vol_rhr)~=length(vol_lhr))
                mess-sprintf('could not read cortical ribbon volume');
                disp(mess)
                i=i+1;
            else
                
                %%  3- Compute Dice coefficient   %%
                %% ------------------------------ %%
                
                n_aseg=0;
                n_ribbon=0;
                n_overlap=0;
                
                sz=size(vol1);
                for u=1:sz(1)
                    for j=1:sz(2)
                        for k=1:sz(3)
                            if (vol1(u,j,k)==42 | vol1(u,j,k)==3) % right & left cerebral cortex in the aseg volume
                                n_aseg = n_aseg+1;
                            end
                            if (vol_lhr(u,j,k)~=0 | vol_rhr(u,j,k)~=0)
                                n_ribbon = n_ribbon+1;
                            end
                            if ((vol_lhr(u,j,k)~=0)  && (vol1(j,u,k)==3))
                                n_overlap = n_overlap+1;
                            end 
                            if ((vol_rhr(u,j,k)~=0)  && (vol1(j,u,k)==42))
                                n_overlap = n_overlap+1;
                            end    
                        end
                    end
                end
                
                d= 2*n_overlap/(n_aseg+n_ribbon);
                D=[D d];
                Isubj=[Isubj i_subj];
                %thres=0.62;  % Thresholds deduced from stats on the training set /space/neo/2/recon/buckner
                pval=compute_pval(d, DD);
                if (pval < th_pval)
                    mess=sprintf('The Dice coefficient is too low for subject %s (D=%3g, pval=%g)', subject_name, d, pval);
                    disp(mess)
                else
                    mess=sprintf('cortical ribbon verification OK for subject %s (D=%3g, pval=%g)', subject_name, d, pval);
                    disp(mess)
                end
                cmdrm=sprintf('rm %s | rm %s', lh_ribbon, rh_ribbon);
                u=unix(cmdrm);
            end
        end
    else
        i=i+1;
    end
end


% subfunction compute_pval() %
function [p_inf]=compute_pval(val, D)
%load('/space/okapi/3/data/laurence/ADF/ribbon/Dice_ribbon.mat'); %loads D
pas=0.05;
x=0:pas:1;
[h] = hist(D,x);
p = h/sum(h);
dinf=find(x<=val);
xinf=x(dinf);
pinf=p(1:length(xinf));
if (val>=0 & length(xinf) >1 )
    p_inf=trapz(xinf,pinf)/pas;
elseif (val>=0 & (length(xinf)<2))
    pas2=pas/10;
    x2=0:pas2:1;
    [h2] = hist(D,x2);
    p2 = h2/sum(h2);
    dinf2=find(x2<=val);
    xinf2=x2(dinf2);
    pinf2=p2(1:length(xinf2));
   if(length(xinf2)>1)
        p_inf=trapz(xinf2,pinf2)/pas2;
    else
        p_inf=0;
    end 
else
    p_inf=0;
end
