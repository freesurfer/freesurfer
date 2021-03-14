function [D, Isubj]=pons_cut_dir_adf(dirname,th_pval)
%         For all the subjects in the directory "dirname": 
%           Computes the Dice coefficients D=2Nab/Na+Nb 
%                           where:
%   Na is the volume of the Cerebellum obtrained trough the volume-based labeling
%     Nb is the volume "filled" obtained from the surface-based stream
%              Nab is the volume of the overlap
%


%
% pons_cut_dir_afd.m
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
    msg=sprintf('USAGE: [D]=cc_cut_dir_adf(dirname, th_pval)');
    disp(msg)
end

%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end

%%% Load stats from the Buckner data set %%%
%stat_file='/space/okapi/3/data/laurence/ADF/cutting_planes/PonsCutDice.adf';
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
DDpons=fscanf(fid, '%g');

files=dir(dirname);

D=[];
Isubj=[];
for i=1:length(files)
    SubjectDir=strcat(dirname,'/',files(i).name);
    %disp(SubjectDir)
    CorDir1=strcat(SubjectDir,'/mri/aseg/');
    CorDir2=strcat(SubjectDir,'/mri/filled/');
    d1=dir(CorDir1);
    d2=dir(CorDir2);
    if (length(d1)<3 | length(d2)<3 |(strcmp(files(i).name, '010128_2105')==1) |(strcmp(files(i).name, '010621_vc7110')==1) | ( length(strfind(files(i).name,'0'))==0 &&(length(strfind(files(i).name,'1'))==0 )))
    %if (length(d1)<3 | length(d2)<3 )
        i=i+1; % go to the next subject in the directory
    else
        %%% load the volumes aseg and wm %%%
    [vol1 M1 t1]=load_cor2(SubjectDir,'aseg');
    [vol2 M2 t2]=load_cor2(SubjectDir,'filled');
    
    if ((t1~=1) | (t2~=1))
        i=i+1;
    else
        %%% Compute the Dice coefficient %%%
        d=compute_dice(vol1,vol2);
        D=[D d];
        Isubj=[Isubj i];
        pval=compute_pval(d, DDpons);
        if (pval<th_pval)
            mess=sprintf('Cerebellum & Brain-stem not succesfully removed for subject: %s (dice=%.3g, pval=%.3g)', files(i).name, d, pval);
            disp(mess)
        else 
            mess=sprintf('Cerebellum & Brain-stem succesfully removed for subject: %s (dice=%.3g, pval=%.3g)', files(i).name, d, pval);
            disp(mess)
        end
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
                if (vol1(i,j,k)==7 | vol1(i,j,k)==46 | vol1(i,j,k)==16)
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
    disp(mess)
end

% subfunction compute_pval() %
function [p_sup]=compute_pval(val, Dpons)
%load('/space/okapi/3/data/laurence/ADF/cutting_planes/Dice_pons.mat'); %loads Dpons
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
