function [Dvol,I]=check_ROI_dir(Dirname, th_pval)
%
%                   For all the subjects in a directory: 
%           check if the size of 20 ROIs is within the normal range
% The 20 following ROIs are checked: Left-Lateral-Ventricle Right-Lateral-Ventricle 
%   Left-Hippocampus Right-Hippocampus Left-Thalamus-Proper Right-Thalamus-Proper 
%     Left-Caudate Right-Caudate Left-Putamen Right-Putamen Left-Pallidum 
% Right-Pallidum Left-Amygdala Right-Amygdala Left-Accumbens-area Right-Accumbens-area 
%      3rd-Ventricle 4th-Ventricle Left-Inf-Lat-Vent Right-Inf-Lat-Vent) 
%


%
% subcortical_labeling_dir_afd.m
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
    msg=sprintf('USAGE: check_ROI_dir(Dirname)');
    disp(msg)
end

ROI_Labels=[4 43 17 53 10 49 11 50 12 51 13 52 18 54 26 58 14 15 5 44];
labmapfile=('/space/lyon/1/fsdev/freesurfer_dev');
labmap=strcat(labmapfile,'/','tkmeditColorsCMA');
[label name val1 val2 val3 val4]=textread(labmap,'%d %s %d %d %d %d',89); 

%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end

% Load stats from the Buckner data set %
%stat_file='/space/okapi/3/data/laurence/ADF/subcortical_labeling/SubcorticalLabelingPercentVol.adf';
stat_file=stract(fsafdDir, '/SubcorticalLabelingPercentVol.adf');
fid=fopen(stat_file);
if(fid==-1)
    mess=sprintf('Could not find %s', stat_file);
    error(mess)
end
while(strfind(fgetl(fid), '#'))
    pos=ftell(fid);
end
fseek(fid, pos, 'bof');
Dtmp=fscanf(fid, '%g');
nrow=length(Dtmp)/20;
DD=(reshape(Dtmp, [20 nrow]))';

%M=[0.3644 0.3111 0.3583 0.3526 0.6776 0.6639 0.3042 0.3070 0.4769 0.4859 0.1399 0.1336 0.1729 0.1739 0.0672 0.0589 0.0532 0.2118 0.0330 0.0341];
%Sd=[0.2080 0.1651 0.0484 0.0395 0.0182 0.0035 0.0058 0.0186 0.0153 0.0104 0.0222 0.0198 0.0113 0.0046 0.0066 0.0001 0.0295 0.0233 0.0173 0.0326];

files=dir(Dirname);
nsup=[ 5 5 0.4  0.4  0.76 0.77 0.52 0.46 0.65  0.56 0.2   0.22 0.2  0.2  0.088 0.08  0.5  0.35  0.35  0.35];
ninf=[ 0 0 0.19 0.14 0.5  0.47 0.2  0.18 0.31  0.3  0.12  0.1  0.06 0.06 0.03  0.025 0.03 0.06  0.01  0.01]; 
load_flag=0;
Dvol=[];
I=[];

for i=5:328
    SubjectDir=strcat(Dirname,'/',files(i).name);
    disp(files(i).name)
    CorDir=strcat(SubjectDir,'/mri/aseg/');
    d=dir(CorDir);
    %if (length(d)<3 | ( length(strfind(files(i).name,'0'))==0))
    if (length(d)<3 )
        aseg_vol=strcat(SubjectDir,'/mri/aseg.mgz');
            if(~exist(aseg_vol))
                i=i+1; % go to the next subject in the directory
            else
                vol=load_mgh(aseg_vol);
                load_flag=1;
            end
    else
        [vol mat]=load_cor2(SubjectDir,'aseg');
        load_flag=1;
    end
    if(load_flag==1)    
        %%% Compute the volumes of the ROIs and the corresponding %%%
        %%% percentages of the total brain volume %%%
        c=0;
        y=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
        sz=size(vol);
        for u=1:sz(1)
            for j=1:sz(2)
                for k=1:sz(3)
                    if (vol(u,j,k)~=0)
                        c=c+1;
                    end
                    for l=1:length(ROI_Labels)
                        if (vol(u,j,k)==ROI_Labels(l))
                            y(l)=y(l)+1;
                        end
                    end
                end
            end
        end
        bv=c;
        for v=1:20
            y(v)=y(v)/bv*100;
        end
        Dvol=[Dvol ; y];
        I=[I i];
        vol=[];
        count=0;
        %disp(SubjectDir)
        for j=3:length(ROI_Labels)
            labelname=char(name(ROI_Labels(j)+1));
            if(y(j)<0.009)
                msg=sprintf('The volume of the %s (%.3g%% of the brain) is abnormally low',labelname,y(j));
                disp(msg)
                count=1;
            else
                [pinf, psup]=compute_pval(j,y, DD);
                if ( pinf<th_pval | psup<th_pval )
                    %msg=sprintf('%s: %s represents %.3g%% of the brain (normalrange: [%.3g .. %.3g] )', files(i).name,labelname,y(j),ninf(j),nsup(j));
                    msg=sprintf('%s: %s represents %.3g%% of the brain (normalrange: [%.3g .. %.3g] )', files(i).name,labelname,y(j),ninf(j),nsup(j));
                    disp(msg)
                    count=1;
                end
            end
        end
        if (count==0)
            msg=sprintf('Subject %s is normal',files(i).name);
            disp(msg)
        end
    end
end

% subfunction compute_pval() %
function [p_inf, p_sup]=compute_pval(llabel,M,D)
%load('/space/okapi/3/data/laurence/ADF/subcortical_labeling/PercentVol_labels.mat'); % loads D
pas=0.01;
x=0:pas:1;
% jackknife ?
[h] = hist(D(:,llabel),x);
p = h/sum(h);
dinf=find(x<=M(llabel));
dsup=find(x>M(llabel));
xinf=x(dinf);
xsup=x(dsup);
pinf=p(1:length(xinf));
psup=p(length(x)-length(xsup)+1:end);
if (M(llabel)>=0 & length(xinf) >1 & length(xsup) >1)
    p_inf=trapz(xinf,pinf)/pas;
    p_sup=trapz(xsup,psup)/pas;
elseif (M(llabel)>=0 & (length(xinf)<2 |  length(xsup)<2))
    pas2=pas/5;
    x2=0:pas2:1;
    [h2] = hist(D(:,llabel),x2);
    p2 = h2/sum(h2);
    dinf2=find(x2<=M(llabel));
    dsup2=find(x2>M(llabel));
    xinf2=x2(dinf2);
    xsup2=x2(dsup2);
    pinf2=p2(1:length(xinf2));
    psup2=p2(length(x2)-length(xsup2)+1:end);
    p_inf=trapz(xinf2,pinf2)/pas2;
    p_sup=trapz(xsup2,psup2)/pas2;
else
    p_inf=0;
    p_sup=0;
end








