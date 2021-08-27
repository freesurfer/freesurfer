function [y]=check_ROI(SubjectDir,th_pval)
%
%   For one subject: check if the size of 20 ROIs is within the normal range
% The 20 following ROIs are checked: Left-Lateral-Ventricle Right-Lateral-Ventricle 
%   Left-Hippocampus Right-Hippocampus Left-Thalamus Right-Thalamus 
%     Left-Caudate Right-Caudate Left-Putamen Right-Putamen Left-Pallidum 
% Right-Pallidum Left-Amygdala Right-Amygdala Left-Accumbens-area Right-Accumbens-area 
%      3rd-Ventricle 4th-Ventricle Left-Inf-Lat-Vent Right-Inf-Lat-Vent) 
%
%     th_pval is the threshold for the p_val 
%     training set: '/space/neo/2/recon/buckner/'
%


%
% subcortical_labeling_afd.m
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
    msg=sprintf('USAGE: [y]=check_ROI(SubjectDir, th_pval)');
    disp(msg)
end

ROI_Labels=[4 43 17 53 10 49 11 50 12 51 13 52 18 54 26 58 14 15 5 44];
labmapfile=('/space/lyon/1/fsdev/freesurfer_dev');
labmap=strcat(labmapfile,'/','tkmeditColorsCMA');
[label name val1 val2 val3 val4]=textread(labmap,'%d %s %d %d %d %d',89); 

%%% load the volume aseg %%%
AsegDir=strcat(SubjectDir,'/mri/aseg');
if (length(dir((AsegDir)))<3)
    aseg_vol=strcat(SubjectDir,'/mri/aseg.mgz');
    if(exist(aseg_vol))
        vol=load_mgh(aseg_vol);
    end
else
    [vol mat]=load_cor2(SubjectDir,'aseg');
end

%%% Compute the volumes of the ROIs and the corresponding %%%
%%% percentages of the total brain volume %%%
c=0;
y=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
sz=size(vol);
for i=1:sz(1)
    for j=1:sz(2)
        for k=1:sz(3)
            if (vol(i,j,k)~=0)
                c=c+1;
            end
            for l=1:length(ROI_Labels)
                if (vol(i,j,k)==ROI_Labels(l))
                    y(l)=y(l)+1;
                end
            end
        end
    end
end
bv=c;
for u=1:20
    y(u)=y(u)/bv*100;
end
vol=[];
count=0;
nsup=[ 5 5 0.4  0.4  0.76 0.77 0.52 0.46 0.65  0.56 0.2   0.22 0.2  0.2  0.088 0.08  0.5  0.35  0.35  0.35];
ninf=[ 0 0 0.19 0.14 0.5  0.47 0.2  0.18 0.31  0.3  0.12  0.1  0.06 0.06 0.03  0.025 0.03 0.06  0.01  0.01]; 
for u=3:length(ROI_Labels)
    labelname=char(name(ROI_Labels(u)+1));
    if(y(u)<0.009)
        msg=sprintf('The volume of the %s (%.3g%% of the brain) is abnormally low',labelname,y(u));
        disp(msg)
        count=count+1;
    else
        [Pinf,Psup]=compute_pval(u,y);
        if ( Pinf<th_pval | Psup < th_pval )
            msg=sprintf('The volume of the %s (%.3g%% of the brain) is out of the normal range [%.3g..%.3g]',labelname,y(u),ninf(u), nsup(u));
            disp(msg)
            count=count+1;
        end 
    end
end
if (count==0)
    msg=sprintf('The volume is normal');
    disp(msg)
end


% subfunction compute_pval() %
function [p_inf, p_sup]=compute_pval(llabel,M)
%load('/space/okapi/3/data/laurence/ADF/subcortical_labeling/PercentVol_labels.mat'); % loads D
% Load stats from the Buckner data set %
%stat_file='/space/okapi/3/data/laurence/ADF/subcortical_labeling/SubcorticalLabelingPercentVol.adf';
%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end
stat_file=strcat(fsafdDir, '/SubcorticalLabelingPercentVol.adf');
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
D=(reshape(Dtmp, [20 nrow]))';
pas=0.01;
x=0:pas:1;
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



