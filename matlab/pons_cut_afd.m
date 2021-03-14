function [D]=pons_cut_adf(subject, th_pval)
% For the subject "subject": computes the Dice coefficient D=2*Nab/(Na+Nb)
%  Na is the volume of the Cerebellum+Brain-stem obtrained trough the volume-based labeling
%     Nb is the volume "filled" obtained from the surface-based stream
%              Nab is the volume of the overlap
%


%
% pons_cut_afd.m
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
    msg=sprintf('USAGE: wm_seg_adf(Subject, th_pval)');
    disp(msg)
end

SubjectDir=strcat(subject,'/');
CorDir1=strcat(subject,'/mri/aseg/');
CorDir2=strcat(subject,'/mri/filled/');
d1=dir(CorDir1);
d2=dir(CorDir2);

if (length(d1)<3)
    mess1=sprintf('Cannot find the coronals of the volume aseg for the subject %s',subject);
    error(mess1)
elseif (length(d2)<3)
    mess2=sprintf('Cannot find the coronals of the volume filled for the subject %s',subject);
    error(mess2)
else
    
    %%% load the aseg and wm volumes %%%
    [vol1 M t]=load_cor2(SubjectDir,'aseg');
    [vol2]=load_cor2(SubjectDir,'filled');
    %%% computes the Dice coefficient %%%
    D=compute_dice(vol1,vol2);
    pval=compute_pval(D);
end

if (pval<th_pval)
    msg=sprintf('Cerebellum & Brain-stem not successfully removed (Dice=%.3g, pval_sup=%.3g)', D, pval);
    disp(msg)
else 
    mess=sprintf('Cerebellum & Brain-stem successfully removed (Dice=%.3g, pval_sup=%.3g)', D, pval);
    disp(mess)
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
                if (vol1(i,j,k)==7 | vol1(i,j,k)==46 | vol1(i,j,k)==16 )
                    Na=Na+1;
                end
                if (vol2(i,j,k)~=0)
                    Nb=Nb+1;
                end
                if ((vol1(i,j,k)==7 | vol1(i,j,k)==46 | vol1(i,j,k)==16) && (vol2(i,j,k)~=0)) 
                    Nab=Nab+1;
                end
            end
        end
    end
    n=2*Nab/(Na+Nb);
else
    mess=sprintf('The wm and aseg volumes do not have the same size');
    n=0;
    error(mess)
end


% subfunction compute_pval() %
function [p_sup]=compute_pval(val)
%load('/space/okapi/3/data/laurence/ADF/cutting_planes/Dice_pons.mat'); %loads Dpons
%Load Buckner stats %
%stat_file='/space/okapi/3/data/laurence/ADF/cutting_planes/PonsCutDice.adf';
%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end
stat_file=strcat(fsafdDir, '/PonsCutDice.adf');
fid=fopen(stat_file);
if(fid==-1)
    mess=sprintf('Could not find %s', stat_file);
    error(mess)
end
while(strfind(fgetl(fid), '#'))
    pos=ftell(fid);
end
fseek(fid, pos, 'bof');
Dpons=fscanf(fid, '%g');
pas=0.05;
x=0:pas:1;
[h] = hist(Dpons,x);
p = h/sum(h);
dsup=find(x>val);
xsup=x(dsup);
psup=p(length(x)-length(xsup)+1:end);
if (val>=0 & length(xsup) >1)
    p_sup=trapz(xsup,psup)/pas;
elseif (val>=0 & (length(xsup)<2))
    pas2=pas/10;
    x2=0:pas2:1;
    [h2] = hist(Dpons,x2);
    p2 = h2/sum(h2);
    dsup2=find(x2>val);
    xsup2=x2(dsup2);
    psup2=p2(length(x2)-length(xsup2)+1:end);
    if(length(xsup2)>1)
        p_sup=trapz(xsup2,psup2)/pas2;
    else
        p_sup=0;
    end  
else
    p_sup=0;
end
