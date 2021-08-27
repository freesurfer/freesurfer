function [probas, Isubj,nf]=taldir(dirname, th_pval, DirTable)
%
% Computes the probability of the Talairach transform matrices
%    of all the subjects found in the directory "dirname".  
%  Uses the mean vector and covariance matrix obtained with talairachin_table.m from 
%     the data set (default data set: /space/neo/2/recon/buckner)
%  Uses th_pval as a threshold for the p-values to detect unlikely transform  matrices
%


%
% talairaching_dir_afd.m
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



if (nargin<2 | nargin>3)
    msg=sprintf('USAGE: [probas, Isubj , nf]=taldir(SubjectDir, th_pval, <DirTable>)');
    disp(msg)
end

files=dir(dirname);
nf=files;
probas=[];
Isubj=[];
Pval=[];
count=0;
probas_flag=1;

if(nargin==3)
    fsafdDir=DirTable;
else
    %%% Get the default tables'directory %%%
    if(getenv('FREESURFER_HOME'))
        fsh=getenv('FREESURFER_HOME');
        fsafdDir=strcat(fsh, '/fsafd');
    else
        error(sprintf('Impossible to find FREESURFER_HOME\n'));
    end
end

%%% 1x9 mean vector obtained from the training set %%%
%load('/space/okapi/3/data/laurence/ADF/talairaching/transfo_param_mean.mat');
%file_mu='/space/okapi/3/data/laurence/ADF/talairaching/TalairachingMean.adf'; %Loads mu
file_mu=strcat(fsafdDir, '/TalairachingMean.adf');
fi=fopen(file_mu);
pos=0;
if(fi==-1)
    mess=sprintf('Could not find %s', file_mu);
    error(mess)
else
    while(strfind(fgetl(fi), '#'))  % skip the header
        pos=ftell(fi);
    end
    fseek(fi, pos, 'bof');
    mu=(fscanf(fi, '%g'))';
    fclose(fi);
end

%%% 9x9 covariance matrix obtained from the training set %%%
%load('/space/okapi/3/data/laurence/ADF/talairaching/transfo_param_regularizedCov2.mat'); %loads sigma
%sigma_file='/space/okapi/3/data/laurence/ADF/talairaching/TalairachingCovariance.adf'; 
sigma_file=strcat(fsafdDir, '/TalairachingCovariance.adf');
fis=fopen(sigma_file);
pos=0;
if(fis==-1)
    mess=sprintf('Could not find %s', sigma_file);
    error(mess)
else
    while(strfind(fgetl(fis), '#'))  % skip the header
        pos=ftell(fis);
    end
    fseek(fis, pos, 'bof');
    sig=fscanf(fis, '%g');
    sigma=reshape(sig, [9,9]);
    fclose(fis);
end

%%% Probabilities of the transform matrices  %%%
%stat_file='/space/okapi/3/data/laurence/ADF/talairaching/TalairachingProbas.adf';
stat_file=strcat(fsafdDir, '/TalairachingProbas.adf');
fid=fopen(stat_file);
pos=0;
if(fid==-1)
    mess=sprintf('Could not find %s... compute probas for the first time', stat_file);
    disp(mess)
    probas_flag=0;
else
    while(strfind(fgetl(fid), '#'))  % skip the header
        pos=ftell(fid);
    end
    fseek(fid, pos, 'bof');
    yy=fscanf(fid, '%g');
    fclose(fid);
end

for i=1:(length(files))
    s=strcat(dirname,'/',files(i).name);
    ttfile=strcat(s,'/mri/transforms/talairach.xfm');
    fid=fopen(ttfile, 'r');
    if ( (fid ~= -1) && ( length(strfind(files(i).name,'0'))>=1 | (length(strfind(files(i).name,'1'))>=1 )))
        while feof(fid) == 0
            linef=fgetl(fid);
            nb=findstr(linef, 'Linear_Transform');
            if nb == 1
                pos=ftell(fid);
                break
            end    
        end
        A=(fscanf(fid, '%g',12))';
        A=[A(1) A(2) A(3) A(5) A(6) A(7) A(9) A(10) A(11)]; % if trans=1, the translation parameters are not taken into account
        status=fclose(fid); 
        p=mvnpdf(A,mu,sigma);
        probas=[probas p];
        Isubj=[Isubj i];
        nf(length(probas)).name=files(i).name;
        if(probas_flag)
            [pinf]=compute_pval(p, yy);
            Pval=[Pval pinf];
            if (pinf < th_pval)
                %mess=sprintf('Talairach Transform: %s failed (%g)', files(i).name, p);
                mess=sprintf('Talairach Transform: %s failed, (p=%g pval=%g)', files(i).name, p, pinf);
                disp(mess)
                count=count+1;
            else
                mess2=sprintf('Talairach Transform: %s OK, (p=%g pval=%g)', files(i).name, p, pinf);
                disp(mess2)
            end
        end
    else 
        i=i+1;
        if (fid ~= -1)
            status=fclose(fid); 
        end
    end
end

if (length(probas) == 0)
    messdir=sprintf('No subject found in this directory');
    disp(messdir)
end
messcount=sprintf('Number of unlikely Talairach transforms : %d', count);
disp(messcount)
% if(~probas_flag)
%     outprobas=strcat(DirTable, '/TalairachingProbas.adf');
%     save(outprobas, 'probas', '-ASCII');
% end
%y=probas;
%save('/space/okapi/3/data/laurence/talairaching/transfo_param_probas.mat', 'y');

% subfunction compute_pval() %
function [p_inf]=compute_pval(val,y)
%load('/space/okapi/3/data/laurence/ADF/talairaching/transfo_param_probas.mat'); %loads y
pas=0.05;
x=0:pas:1;
[h] = hist(y,x);
p = h/sum(h);
dinf=find(x<=val);
xinf=x(dinf);
pinf=p(1:length(xinf));
if (val>=0 & length(xinf) >1 )
    p_inf=trapz(xinf,pinf)/pas;
elseif (val>=0 & (length(xinf)<2))
    pas2=pas/10;
    x2=0:pas2:1;
    [h2] = hist(y,x2);
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

