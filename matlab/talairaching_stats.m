function [D,mu,sigma]=talairach_stats_correct(dirname, outdir)
%
% Computes the mean and covariance matrix from a training set
% 
% By default, the 3 translation parameters are not considered 
%     -> mu is a 1x9 vector and sigma a 9x9 matrix


%
% talairaching_stats.m
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



dirname1='/space/neo/2/recon/buckner/';
dirname2='/space/fiat/1/users/buckner/';
dirname3='/space/brainiac/1/users/xhan/freesurfer/Pfizer/LDA/Oct27/';
dirname4='/space/gam/2/users/jjwisco/GSK/cntls/';

mes=sprintf('Talairach stats...');
disp(mes)

if (nargin>0)
    disp(dirname)
    D1=read_talmat(dirname);
else
    D1=read_talmat(dirname1); 
end

D=D1;
D=[D(:,1) D(:,2) D(:,3) D(:,5) D(:,6) D(:,7) D(:,9) D(:,10) D(:,11)];


mu=mean(D);
sigma=cov(D); 

%%% Regularisation of the covariance matrix %%%
[u s v]=svd(sigma);
ds=diag(s);
ds(1:end)=ds(1:end)+0.15;
sigma=u*(diag(ds))*v';
%save('/space/okapi/3/data/laurence/talairaching/transfo_param.mat', 'D');
%save('/space/okapi/3/data/laurence/talairaching/transfo_param_mean.mat', 'mu');
%save('/space/okapi/3/data/laurence/talairaching/transfo_param_regularizedCov.mat', 'sigma');
outmean=strcat(outdir, '/TalairachingMean_tmp.adf');
outcov=strcat(outdir, '/TalairachingCovariance_tmp.adf');
save(outmean, 'mu', '-ASCII');
save(outcov, 'sigma', '-ASCII');

%%%   Subfunction read_talmat(dirname,opt)   %%%
%%% Collect the talairach parameters of all  %%%
%%%   subjects in the directory "dirname"    %%%
function D=read_talmat(dname)
D=0;
files=dir(dname);
for i=1:(length(files))
    s=strcat(dname,'/',files(i).name);
    ttfile=strcat(s,'/mri/transforms/talairach.xfm');
    fid=fopen(ttfile, 'r');
   if ( (fid ~= -1) && ( length(strfind(files(i).name,'0'))>=1 ||(length(strfind(files(i).name,'1'))>=1 )))
       while feof(fid) == 0
           linef=fgetl(fid);
           nb=findstr(linef, 'Linear_Transform');
           if nb == 1
               pos=ftell(fid);
               break
           end    
       end
       A=(fscanf(fid, '%g',12));
       status=fclose(fid);
       if D==0
           D=A';
       else
           D=[D ; A'];
       end
       
   else
        i=i+1;
    end
end
