function [gca] = load_gca(fname)
%
% [gca] = load_gca(fname)
% reads an array of gaussian classifiers
%


%
% load_gca.m
%
% Original Author: Bruce Fischl
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

GCA_VERSION=4.0;
GCA_NO_MRF=1;
GIBBS_NEIGHBORHOOD=6 ;
GIBBS_NEIGHBORS=GIBBS_NEIGHBORHOOD;
MAX_LABELS=4;

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
	 str = sprintf('could not open gca file %s.', fname) ;
	 error(str) ;
end
version = fread(fid, 1, 'float32') ;
if (version < 1 | version > GCA_VERSION)
	 fclose(fid) ;
	 error(sprintf('version %d in file %s incorrect - not a known gca file',version,fname));
end

prior_spacing = fread(fid, 1, 'float32') ;
node_spacing = fread(fid, 1, 'float32') ;

prior_width = fread(fid, 1, 'int32') ;
prior_height = fread(fid, 1, 'int32') ;
prior_depth = fread(fid, 1, 'int32') ;

node_width = fread(fid, 1, 'int32') ;
node_height = fread(fid, 1, 'int32') ;
node_depth = fread(fid, 1, 'int32') ;

ninputs = fread(fid, 1, 'int32') ;
flags = fread(fid, 1, 'int32') ;

disp(sprintf('reading gca file %s (%dx%dx%d), spacing=%d, version %.1f', fname, prior_width,prior_height,prior_depth,prior_spacing,version)) ;

gca = zeros(prior_width*prior_height*prior_depth, 2*MAX_LABELS+1) ;

index = 1 ;
for x=1:node_width
%		disp(sprintf('reading slice %d of %d', x, node_width)) ;
		for y=1:node_height
				for z=1:node_depth
						nlabels = fread(fid, 1, 'int32') ;
						total_training = fread(fid, 1, 'int32') ;
%						gca(index,1) = nlabels ;

						for n=1:nlabels
								label = fread(fid, 1, 'uchar') ;
								mean = fread(fid, 1, 'float32') ;
								var = fread(fid, 1, 'float32') ;
								if (bitand(flags, GCA_NO_MRF))
									continue ;
							  end
								for i=1:GIBBS_NEIGHBORS
										gibbs_nlabels = fread(fid, 1, 'uint32') ;
										for j=1:gibbs_nlabels
												gibbs_label = fread(fid, 1, 'uint32') ;
												gibbs_prior = fread(fid, 1, 'float32') ;
										end
								end
						end
						index = index+1 ;
				end
		end
end


index = 1 ;
for x=1:prior_width
		disp(sprintf('reading slice %d of %d', x, prior_width)) ;
		for y=1:prior_height
				for z=1:prior_depth
						nlabels = fread(fid, 1, 'int32') ;
						total_training = fread(fid, 1, 'int32') ;
						gca(index,1) = nlabels ;

						if (x == prior_width/2 & y==prior_height/2 & z==prior_depth/2)
%							 keyboard ;
						end
						for n=1:nlabels
								label = fread(fid, 1, 'uchar') ;
								prior = fread(fid, 1, 'float32') ;
								if (n <= MAX_LABELS)
									gca(index,2*n) = label ;
									gca(index,2*n+1) = prior ;
								end
						end
						index = index+1 ;
				end
		end
end

fclose(fid) ;

