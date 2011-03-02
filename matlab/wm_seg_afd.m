function [Dice,pval]=wm_seg_adf(subject, th_pval)
% For the subject "subject": computes the Dice coefficient D=2Nab/Na+Nb 
%  where:
%   Na is the volume of the WM obtrained trough the volume-based labeling
%   Nb is the volume of the WM segmented in the surface-based stream
%   Nab is the volume of the overlap
%  Uses the p values
%


%
% wm_seg_afd.m
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


if (nargin<2 | nargin>2)
    msg=sprintf('USAGE: [Dice,pval]=wm_seg_adf(Subject, th_pval)');
    disp(msg)
end


SubjectDir=strcat(subject,'/');
CorDir1=strcat(subject,'/mri/aseg/');
CorDir2=strcat(subject,'/mri/wm/');
d1=dir(CorDir1);
d2=dir(CorDir2);

if (length(d1)<3)
    mess1=sprintf('Cannot find the coronals of the volume aseg for the subject %s',subject);
    disp(mess1)
elseif (length(d2)<3)
    mess2=sprintf('Cannot find the coronals of the volume wm for the subject %s',subject);
        disp(mess2)
else
    
    %%% load the aseg and wm volumes %%%
    [vol1 M t]=load_cor2(SubjectDir,'aseg');
    [vol2]=load_cor2(SubjectDir,'wm');
    
    %%% computes the Dice coefficient %%%
    Dice=compute_dice(vol1,vol2);
    pval=compute_pval(Dice);
    if (pval<th_pval) % (th~=D=0.72)threshold found from the training set '/space/neo/2/recon/buckner'
        msg=sprintf('The dice coefficient is too low (D:%.3g, pval: %.3g)', Dice, pval);
        disp(msg)
    else
        mess=sprintf('WM segmentation OK (dice coefficient = %.3g, pval: %.3g)',Dice, pval);
        disp(mess)
    end
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
                if (vol1(i,j,k)==2 | vol1(i,j,k)==41 ) % aseg volume
                    Na=Na+1;
                end
                if (vol2(i,j,k)~=0) % wm volume
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
function [p_inf]=compute_pval(val)
%load('/space/okapi/3/data/laurence/ADF/wm_segmentation/WM_dice_buckner.mat'); %loads D
%stat_file='/space/okapi/3/data/laurence/ADF/wm_segmentation/WMsegmentationDice.adf';
%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end
stat_file=strcat(fsafdDir, '/WMsegmentationDice.adf');
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

