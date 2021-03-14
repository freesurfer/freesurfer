function [Dice, Pval, Isubj]=wm_seg_dir_adf(dirname, th_pval)
%         For all the subjects in the directory "dirname": 
%           Computes the Dice coefficients D=2Nab/Na+Nb 
%                           where:
%   Na is the volume of the WM obtained trough the volume-based labeling
%     Nb is the volume of the WM segmented in the surface-based stream
%              Nab is the volume of the overlap
%
%  Uses the p values (NB, Dices coefficients take some time to compute...)
%


%
% wm_seg_dir_afd.m
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



if (nargin<2 | nargin>2)
    msg=sprintf('USAGE: [Dice,Pval, Isubj]=wm_seg_dir_adf(dirname, th_pval)');
    disp(msg)
end

%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end

% Load the Buckner stats % 
%stat_file='/space/okapi/3/data/laurence/ADF/wm_segmentation/WMsegmentationDice.adf';
stat_file=strcat(fsafdDir, '/WMsegmentationDice.adf');
fid=fopen(stat_file);
if(fid==-1)
    mess=sprintf('Could not find %s', stat_file);
    error(mess)
end
while(strfind(fgetl(fid), '#'))  % skip the header
    pos=ftell(fid);
end
fseek(fid, pos, 'bof');
DD=fscanf(fid, '%g');


files=dir(dirname);

Dice=[];
P=[];
Isubj=[];

for i=1:length(files)
    SubjectDir=strcat(dirname,'/',files(i).name);
    CorDir1=strcat(SubjectDir,'/mri/aseg/');
    CorDir2=strcat(SubjectDir,'/mri/wm/');
    d1=dir(CorDir1);
    d2=dir(CorDir2);
    if (length(d1)<3 | length(d2)<3 | (strcmp(files(i).name,'010611_vc7044')==1) | ( length(strfind(files(i).name,'0'))==0 &&(length(strfind(files(i).name,'1'))==0 )))
    %if (length(d1)<3 | length(d2)<3 )
        i=i+1; % go to the next subject in the directory
    else
        %%% load the volumes aseg and wm %%%
        [vol1 M1 t1]=load_cor2(SubjectDir,'aseg');
        [vol2 M2 t2]=load_cor2(SubjectDir,'wm');
        
        if ((t1~=1) | (t2~=1))
            i=i+1;
        else
            %%% Compute the Dice coefficient %%%
            d=compute_dice(vol1,vol2);
            pval=compute_pval(d, DD);
            %disp(sprintf('subj %s, d=%.3g, pval=%.3g', files(i).name, d, pval));
            if (pval<th_pval)
                msg=sprintf('Wrong wm segmentation for subject %s, (d=%.3g, pval=%.3g)', files(i).name, d, pval);
                disp(msg)
            else
                msg2=sprintf('wm segmentation OK for subject %s, (d=%.3g, pval=%.3g)', files(i).name, d, pval);
                disp(msg2)
            end
            
            Dice=[Dice d];
            P=[P pval];
            Isubj=[Isubj i];
            
        end
    end
    vol1=[];
    vol2=[];
end

%%% sub-function [n]=compute_dice(vol1,vol2) %%%

function [n]=compute_dice(vol1,vol2)
 
Na=0;
Nb=0;
Nab=0;
sz=size(vol1);
sw=size(vol2);
if (sz==sw)
    for i=1:sz(1)
        for j=1:sz(2)
            for k=1:sz(3)
                if (vol1(i,j,k)==2 | vol1(i,j,k)==41) % Does not include the cerebellum WM ( keep this ???)
                    Na=Na+1;
                end
                if (vol2(i,j,k)~=0)
                    Nb=Nb+1;
                end
                if ((vol1(i,j,k)==2 | vol1(i,j,k)==41 ) && (vol2(i,j,k)~=0)) % Does not include the cerebellum WM
                    Nab=Nab+1;
                end
            end
        end
    end
    n=2*Nab/(Na+Nb);
else
    mess=sprintf('The wm and aseg volumes do not have the same size');
    n=0;
    disp(mess)
end

% subfunction compute_pval() %
function [p_inf]=compute_pval(val, D)
%load('/space/okapi/3/data/laurence/ADF/wm_segmentation/WM_dice_buckner.mat'); %loads D
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
