function [Dice]=ribbon_adf(subject, th_pval)

%           For each subject "subject":
% Computes the Dice coefficients measuring the overlap of the
% Cortical Ribbon volume computed 
%         1- from the subcortical labeling   
%         2- as the space between the white and the pial surface 
%
%  Uses the pvalues, try to use rh.ribbon.mgh/lh.ribbon.mgh in the subject/surf directory
%  create these volumes if they don't exist
%


%
% ribbon_afd.m
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




if (nargin<1 | nargin>2)
    msg=sprintf('USAGE: [D]=ribbon_adf(Subject, th_pval)');
    disp(msg)
end

%%  1- Load the aseg volume    %%
%% --------------------------  %%

SubjectDir=strcat(subject,'/');
CorDir1=strcat(subject,'/mri/aseg/');
d1=dir(CorDir1);

if (length(d1)<3)
    mess=sprintf('Cannot find the coronals of the volume aseg for the subject %s',subject);
    error(mess)
else
    [vol1 M]=load_cor2(SubjectDir,'aseg');
    
end

%%  2- Check if rh.ribbon.mgh / lh.ribbon.mgh exist     %%
%%  in the /surf directory, or create them              %%
%% ---------------------------------------------------  %%

create_ribbon=0;
SurfDir=strcat(SubjectDir, 'surf');
surf=dir(SurfDir);
if(length(surf)<=2)
    mess=sprintf('Directory %s is empty', SurfDir);
    error(mess)
end
if(length(surf)>2)
    if (( ~(exist(strcat(SurfDir, '/rh.ribbon.mgh'))) | ~(exist(strcat(SurfDir, '/lh.ribbon.mgh'))) ) & ( ~(exist(strcat(SurfDir, '/rh.ribbon.mgz'))) | ~(exist(strcat(SurfDir, '/lh.ribbon.mgz'))) ) )
        mess=sprintf('Cannot find the cortical ribbon for subject %s... try to create it', subject);
        disp(mess)
        create_ribbon=1;
    end
end
if(length(surf)>2 & create_ribbon==1)
    if (~(exist(strcat(SurfDir, '/rh.white'))) | ~(exist(strcat(SurfDir, '/lh.white'))) | ~(exist(strcat(SurfDir, '/rh.thickness'))) | ~(exist(strcat(SurfDir, '/lh.thickness'))))
        mess=sprintf('Cannot find the white surfaces and thicknesses for subject %s', subject);
        error(mess)
    end
end

if(create_ribbon==1)
    % b- Generate a volume registration file 
    a=strread(subject, '%s', 'delimiter', '/');
    subject_name=char(a(length(a)));
    %disp(subject_name)
    reg_file=sprintf('~/reg.dat');
    fid=fopen(reg_file, 'w');
    fprintf(fid, '%s\n', subject_name);
    fprintf(fid, '1\n1\n1\n1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1\n%s\n', 'round');
    fclose(fid);
    
    
    % c- Fill ribbon in the left & right hemispheres
    
    template_vol=strcat(subject, '/mri/orig/');
    lh_ribbon=sprintf('~/ribbon_lh.mgz'); % temporary location for the created ribbon volumes
    rh_ribbon=sprintf('~/ribbon_rh.mgz');
    %lh_ribbon=strcat(subject,'/surf/lh.ribbon.mgh'); % definitiv location for the created ribbon volumes
    %rh_ribbon=strcat(subject,'/surf/rh.ribbon.mgh');
    
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
            error(msg);
        end
    end
    
    cmd1=sprintf('~/work/dev/mri_surf2vol/mri_surf2vol --mkmask --hemi lh --fillribbon --template %s --outvol %s --volreg %s --sd %s', template_vol, lh_ribbon, reg_file, SubjectsDir);
    cmd2=sprintf('~/work/dev/mri_surf2vol/mri_surf2vol --mkmask --hemi rh --fillribbon --template %s --outvol %s --volreg %s --sd %s', template_vol, rh_ribbon, reg_file, SubjectsDir);
    
    c1=unix(cmd1);
    c2=unix(cmd2);
else
    if(exist(strcat(SurfDir, '/lh.ribbon.mgh')) & exist(strcat(SurfDir, '/rh.ribbon.mgh')) )
        lh_ribbon=strcat(SurfDir, '/lh.ribbon.mgh');
        rh_ribbon=sprintf(SurfDir, '/rh.ribbon.mgh');
    end
    if(exist(strcat(SurfDir, '/lh.ribbon.mgz')) & exist(strcat(SurfDir, '/rh.ribbon.mgz')) )
        lh_ribbon=strcat(SurfDir, '/lh.ribbon.mgz');
        rh_ribbon=sprintf(SurfDir, '/rh.ribbon.mgz');
    end
end

% d- Load the right & left cortical ribbon

if( (exist(lh_ribbon)==0) | (exist(rh_ribbon)==0) )
    mess=sprintf('Could not find cortical ribbon volume');
    error(mess)
else
    [vol_lhr]=load_mgh(lh_ribbon);
    [vol_rhr]=load_mgh(rh_ribbon);

    %%  3- Compute Dice coefficient   %%
    %% ------------------------------ %%

    n_aseg=0;
    n_ribbon=0;
    n_overlap=0;

    if ( (size(vol1)) ~= (size(vol_lhr)) | (size(vol1)) ~= (size(vol_rhr)) )
        msg=sprintf('The different volumes do not have the same size');
        error(msg);
    else
        sz=size(vol1);
        for i=1:sz(1)
            for j=1:sz(2)
                for k=1:sz(3)
                    if (vol1(i,j,k)==42 | vol1(i,j,k)==3) % right & left cerebral cortex in the aseg volume
                        n_aseg = n_aseg+1;
                    end
                    if (vol_lhr(i,j,k)~=0 | vol_rhr(i,j,k)~=0)
                        n_ribbon = n_ribbon+1;
                    end
                    if ((vol_lhr(i,j,k)~=0)  && (vol1(j,i,k)==3))
                        n_overlap = n_overlap+1;
                    end 
                    if ((vol_rhr(i,j,k)~=0)  && (vol1(j,i,k)==42))
                        n_overlap = n_overlap+1;
                    end    
                end
            end
        end
    end

    Dice = 2*n_overlap/(n_aseg+n_ribbon);
    pval=compute_pval(Dice);
    %thres=0.0092;  % P(Dice=0.6) for the Buckner training set
    %thres=0.62 %for the Dice value...
    if (pval < th_pval)
        mess=sprintf('The Dice coefficient is too low (dice:%3g, pval:%3g)', Dice, pval);
        disp(mess)
    else
        mess=sprintf('cortical ribbon verification OK (dice:%3g, pval:%3g)', Dice, pval);
        disp(mess)
    end
    cmdrm=sprintf('rm %s | rm %s', lh_ribbon, rh_ribbon);
    %u=unix(cmdrm);
end                
    

% subfunction compute_pval() %
function [p_inf]=compute_pval(val)
%load('/space/okapi/3/data/laurence/ADF/ribbon/Dice_ribbon.mat'); %loads D
% Load stats from the Buckner data set %
%stat_file='/space/okapi/3/data/laurence/ADF/ribbon/RibbonPlacementDice.adf';
%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end
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
D=fscanf(fid, '%g');
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



    
    

        
    
  
   
