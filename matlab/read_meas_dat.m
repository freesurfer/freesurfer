function varargout = read_meas_dat(filename, options)
%READ_MEAS_DAT  read in VB13A-style "meas.dat"
%
% data = read_meas_dat(filename,<options>)
%
% [data, phascor1d, phascor2d, noise] = read_meas_dat(filename)
%
% [data, phascor1d, phascor2d, noise, patrefscan, patrefscan_phascor] = read_meas_dat(filename)
%
% options is a structure which allows changing of some of the defaults:
%   options.ReverseLines - set to 0 or 1 (default is 1)
%   options.ApplyFFTScaleFactors - set to 0 or 1 (default is 0)
%   options.CanonicalReorderCoilChannels - set to 0 or 1 (default 1)
% If a field does not exist, then the default is used.
% options is optional :).


% (based on the FAMOUS "read_mdh_adc.m" by anders, then mukund, then andre.)

% 2006/dec/04: added support for EPI data
%                - (PHASCOR) extract phase correction lines
%                - (REFLECT) reorder samples from reflected lines

% 2006/dec/06: added support for iPAT data
%                - (NOISEADJSCAN)
%                - (PATREFSCAN)

% 2007/jan/01: added support for 3D EPI phase correction lines

% 2007/mar/29: added canonical preamp-based coil ordering (cf. ICE code)

% 2007/mar/30: implemented Anastasia's suggestions
%                - (PATREFANDIMASCAN) for non-EPI iPAT data

% 2007/may/02: added support for phase stabilization scans
%                - (REFPHASESTABSCAN)
%                - (PHASESTABSCAN)

% jonathan polimeni <jonnyreb@padkeemao.nmr.mgh.harvard.edu>, 10/04/2006
% $Id: read_meas_dat.m,v 1.3 2007/05/09 17:44:33 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.3 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %------------------------------------------------------------------------%
  % error checking

  if ( ~exist(filename, 'file') ),
    error('file [%s] does not exist', filename);
  end;


  %------------------------------------------------------------------------%
  % the pre-release version installed on Bay4 in 2005 differs slightly from
  % the full-release version, and files must be read in differently
  % depending on the version. yuck!
  IS__VB13_PRE_RELEASE_VERSION = 0;
  IS__VB13A_RELEASE_VERSION = 1;

  DO__MDH_SAVE = 0;

  DO__REVERSE_LINES = 1;
  DO__CANONICAL_REORDER_COIL_CHANNELS = 0;
  DO__APPLY_FFT_SCALEFACTORS = 0;

  % If the options variable is passsed, then override
  if ( exist('options','var') ) 
    if(isfield(options,'ReverseLines'))
      DO__REVERSE_LINES = options.ReverseLines;
      fprintf('DO__REVERSE_LINES %d\n',DO__REVERSE_LINES);
    end
    if(isfield(options,'CanonicalReorderCoilChannels'))
      DO__CANONICAL_REORDER_COIL_CHANNELS = options.CanonicalReorderCoilChannels;
      fprintf('DO__CANONICAL_REORDER_COIL_CHANNELS %d\n',DO__CANONICAL_REORDER_COIL_CHANNELS);
    end
    if(isfield(options,'ApplyFFTScaleFactors'))
      DO__APPLY_FFT_SCALEFACTORS = options.ApplyFFTScaleFactors;
      fprintf('DO__APPLY_FFT_SCALEFACTORS %d\n',DO__APPLY_FFT_SCALEFACTORS);
    end
  end
  
  
  t0 = clock;

  [fp, errstr] = fopen(filename, 'r', 'l');
  if ( fp == -1 ),
    error(errstr);
  end;


  % constants defined in <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h>
  MDH_NUMBEROFEVALINFOMASK   = 2;
  MDH_NUMBEROFICEPROGRAMPARA = 4;

  MDH_FREEHDRPARA = 4;

  % from "MrServers/MrVista/include/Ice/IceDefs.h":
  ICE_RAWDATA_SCALE       = 131072.0;  % 64 ^ 3 / 2
  K_ICE_AMPL_SCALE_FACTOR = 80 * 20 * ICE_RAWDATA_SCALE / 65536;


  % determine size (in bytes) of ascii header files stored in the 'meas.dat'
  % format (i.e., "Config_.evp", "Dicom_.evp", etc) to skip over them all.
  % [note: this works for VB11A also---the first integer is 32, which is
  % the number of bytes to skip at the beginning!]
  data_start = fread(fp, 1, 'uint32');

  % read header into one string for parsing
  header = fscanf(fp, '%c', data_start-4);

  param_list = {'NColMeas', 'NLinMeas', 'NChaMeas', 'NSetMeas', 'NEcoMeas', ...
                'NPhsMeas', 'NRepMeas', 'NSegMeas', 'NParMeas', 'NSlcMeas', ...
                'NIdaMeas', 'NIdbMeas', 'NIdcMeas', 'NIddMeas', 'NIdeMeas', ...
                'NAveMeas'};

  dimensions = cell2struct(cell(length(param_list),1), param_list, 1);
  dim = [];

  % scan through header for each of the ICE dimension values
  for ind = 1:length(param_list),
    param = param_list{ind};

    %%% SPECIAL CASE: "NSegMeas"

    % the number of segments is listed in two places in the header with the
    % field names "NSegMeas" and "NSeg", and for some reason only the "NSeg"
    % field gives the correct number of segments; SO for this field we break
    % with the convention [TODO: fix the bug in the martinos EPI sequence
    % that causes this hack]
    % UPDATE: it appears that "NSegMeas" corresponds to the cumulative
    % number of distinct segments appearing in the loop counters, whereas
    % "NSeg" is the true number of segments that is the same as the
    % number of shots, i.e., they should differ by a factor of 2 when the
    % "OnlineTSE" functor is in use.
    if ( strcmp(param, 'NSegMeas') ),
      param = 'NSeg';
    end;

    % exploit MATLAB regexp machinery to pull out parameter/value pairs
    match = regexp(header, ['(?<param>' param, ').{0,5}\{\s*(?<value>\d*)\s*\}'], 'names');

    % check if no match is found
    if ( isempty(match) ),
      if ( IS__VB13A_RELEASE_VERSION ),
        IS__VB13A_RELEASE_VERSION = 0;
        warning('SIEMENS:IO:versioning', 'missing header data---check ICE version');
      end;
      continue;
    end;

    % consider only last match (there can be as many as three in Config_.evp)
    match = match(end);

    % empty means number of elements in this dimension = 1
    if ( isempty(match.value) ),
      match.value = '1';
    end;

    % save out struct and numerical array
    dim(ind) = str2double(match.value);
    dimensions.(param_list{ind}) = dim(ind);

  end;

  % VB11A workaround hack
  if ( ~IS__VB13A_RELEASE_VERSION ),
    dimensions.NColMeas = 0;
    dimensions.NAveMeas = 1;
    dimensions.NSegMeas = 2;
  end;


  %------------------------------------------------------------------------%
  % extract FFT scalefactors (if utility file is in path)

  if ( DO__APPLY_FFT_SCALEFACTORS && ...
       exist('read_meas_dat__fft_scalefactors', 'file') ),
    fft_scale = read_meas_dat__fft_scalefactors(header);
    disp('applying FFT scale factors.');
  else,
    fft_scale = ones(dimensions.NChaMeas, 1);
    disp('ignoring FFT scale factors.');
  end;


  %------------------------------------------------------------------------%
  % compute canonical coil ordering from coil element strings--- the channel
  % number assigned to each coil is determined by the order in which the
  % channels are selected in SYNGO, and NOT by the preamp channel numbers,
  % so to enforce a consistent ordering across scans that is independent of
  % the operator's coil selections we can sort by the fixed preamp channel
  % strings.

  if ( DO__CANONICAL_REORDER_COIL_CHANNELS ),

    % prepend "asCoilSelectMeas" to avoid confusion with "sCoilSelectUI" list.

    match = regexp(header, ['asCoilSelectMeas\S*\.sCoilElementID\.tCoilID\s*=\s*"(?<string>\w*)"'], 'names', 'once');
    coil_id = match.string;

    match = regexp(header, ['asCoilSelectMeas\S*\.sCoilElementID\.tElement\s*=\s*"(?<string>\w*)"'], 'names');
    coil_element_string = {match(1:end/2).string};

    match = regexp(header, ['asCoilSelectMeas\S*\.lRxChannelConnected\s*=\s*(?<value>\w*)'], 'names');
    coil_selected = str2num(str2mat({match(1:end/2).value}.'));

    [sorted, coil_index] = sort(coil_element_string);
    [sorted, coil_order] = sort(coil_index);

    % to map using index:
    %   dataval(coil_index);
    % to map using order:
    %   coil_order(dataind);

    disp('reordering coil channels to canonical preamp-based ordering.');

    % don't forget to reorder FFT scalefactors!
    fft_scale = fft_scale(coil_index);

  end;


  %------------------------------------------------------------------------%

  % jump past header to beginning of binary data
  fseek(fp, data_start, 'bof');


  % save current file position, then jump to end to calculate remaining data size
  fpos = ftell(fp);
  fseek(fp, 0, 'eof');
  eof = ftell(fp);
  databytes = eof;

  % return to saved file position
  fseek(fp, fpos, 'bof');

  disp(sprintf('binary data after header: %10.2f MB', databytes/2^20));


  % initializations...
  meas_num = 0;
  scan_num = 0;
  echo_num = 0;

  if ( DO__MDH_SAVE ),
    mdh = repmat(read_meas__mdh_struct, [1, total_scans]);
  else,
    mdh = read_meas__mdh_struct;
  end;

  % number of measurements along all 16 dimensions can be found in the
  % 'Config_.evp' file contained within "meas.dat" header; allocate memory
  % for a single measurement data, and grow data matrix with each line (this
  % method is, surprisingly, faster than pre-allocation!)
%  data                    = single(complex(zeros(dimensions.NColMeas, 1)));
  data                    = single([]);

  FLAG_data_reflect       = 0;

  FLAG_phascor1d          = 0;  % indicating presence of any type of phascor line

  data_phascor1d          = single([]);
  FLAG_data_phascor1d     = 0;
  data_phascor2d          = single([]);
  FLAG_data_phascor2d     = 0;
  data_fieldmap           = single([]);
  FLAG_data_fieldmap      = 0;
  noiseadjscan            = single([]);
  FLAG_noiseadjscan       = 0;
  patrefscan              = single([]);
  FLAG_patrefscan         = 0;
  patrefscan_phascor      = single([]);
  FLAG_patrefscan_phascor = 0;
  FLAG_patrefandimascan   = 0;
  phasestabscan           = single([]);
  refphasestabscan        = single([]);
  FLAG_phasestabscan      = 0;

  FLAG_multiple_repetitions = 0;

  FLAG_data_phascor1d_orderswap = 0;


  % convert binary mask into strings
  EvalInfoMask = read_meas__evaluation_info_mask_definition;
  EvalInfoNames = fieldnames(EvalInfoMask);


  % channel index offset (needed only if IS__VB13_PRE_RELEASE_VERSION)
  channel_offset = NaN;

  ACQEND = 0;

  try,
  %------------------------------------------------------------------------%
  % loop through file, peeling off one MDH and ADC line at a time
  while ( ~ACQEND ),

    meas_num = meas_num + 1;

    if ( DO__MDH_SAVE ),
      idx = meas_num;
    else,
      % overwrite mdh(1) with every measurement line
      idx = 1;
    end;

    mdh(idx).ulFlagsAndDMALength        = uint32(fread(fp, 1, 'uint32'));
    mdh(idx).lMeasUID                   = int32(fread(fp, 1, 'int32'));
    mdh(idx).ulScanCounter              = uint32(fread(fp, 1, 'uint32'));
    if ( ~isempty(mdh(idx).ulScanCounter) ), scan_num = mdh(idx).ulScanCounter; end;
    mdh(idx).ulTimeStamp                = uint32(fread(fp, 1, 'uint32')*2.5); % milliseconds
    mdh(idx).ulPMUTimeStamp             = uint32(fread(fp, 1, 'uint32')*2.5); % milliseconds
    mdh(idx).aulEvalInfoMask(1:MDH_NUMBEROFEVALINFOMASK) = uint32(fread(fp, MDH_NUMBEROFEVALINFOMASK, 'uint32'));

    % build 64-bit mask from two 32-bit integers
    mask = uint64(double(...
        bitshift(uint64(mdh(idx).aulEvalInfoMask(2)), 32)) ...
                  + double(mdh(idx).aulEvalInfoMask(1)));

    mdh(idx).sEvalInfoMask = EvalInfoNames(bitget(mask, 1:length(EvalInfoNames))==1);

    mdh(idx).ushSamplesInScan           = uint16(fread(fp, 1, 'uint16'));
    mdh(idx).ushUsedChannels            = uint16(fread(fp, 1, 'uint16'));

    if (1),
      %  sLoopCounter
      mdh(idx).ushLine                    = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushAcquisition             = uint16(fread(fp, 1, 'uint16'));  % note: acquisition is same as average
      mdh(idx).ushSlice                   = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushPartition               = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushEcho                    = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushPhase                   = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushRepetition              = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushSet                     = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushSeg                     = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushIda                     = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushIdb                     = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushIdc                     = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushIdd                     = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushIde                     = uint16(fread(fp, 1, 'uint16'));
    end;

    if (1),
      % sCutOffData
      mdh(idx).ushPre                     = uint16(fread(fp, 1, 'uint16'));
      mdh(idx).ushPost                    = uint16(fread(fp, 1, 'uint16'));
    end;

    mdh(idx).ushKSpaceCentreColumn      = uint16(fread(fp, 1, 'uint16'));
    mdh(idx).ushDummy                   = uint16(fread(fp, 1, 'uint16'));
    mdh(idx).fReadOutOffcentre          = single(fread(fp, 1, 'float32'));
    mdh(idx).ulTimeSinceLastRF          = uint32(fread(fp, 1, 'uint32'));
    mdh(idx).ushKSpaceCentreLineNo      = uint16(fread(fp, 1, 'uint16'));
    mdh(idx).ushKSpaceCentrePartitionNo = uint16(fread(fp, 1, 'uint16'));

    mdh(idx).aushIceProgramPara(1:MDH_NUMBEROFICEPROGRAMPARA) = uint16(fread(fp, MDH_NUMBEROFICEPROGRAMPARA, 'uint16'));
    mdh(idx).aushFreePara(1:MDH_FREEHDRPARA)                  = uint16(fread(fp, MDH_FREEHDRPARA, 'uint16'));

    if (1),
      % sSliceData
      if (1),
        % sVector
        mdh(idx).flSag            = single(fread(fp, 1, 'float32'));
        mdh(idx).flCor            = single(fread(fp, 1, 'float32'));
        mdh(idx).flTra            = single(fread(fp, 1, 'float32'));
      end;
      mdh(idx).aflQuaternion(1:4) = single(fread(fp, 4, 'float32'));
    end;

    mdh(idx).ushChannelId  = uint16(fread(fp, 1, 'uint16'));
    mdh(idx).ushPTABPosNeg = uint16(fread(fp, 1, 'uint16'));

    if ( IS__VB13_PRE_RELEASE_VERSION ),

      % prior to the full release, a "known bug" in the channel ID numbering
      % conspired to index the channels sequentially BUT began the indexing with
      % a seemingly *random* integer. (twitzel implemented a similar workaround in
      % his code.)

      FIRSTSCANINSLICE = read_meas__extract_flag(mdh(idx).aulEvalInfoMask(1), 'FIRSTSCANINSLICE');
      if ( FIRSTSCANINSLICE && (mdh(idx).ulScanCounter == 1) ),

        if ( isnan(channel_offset) ),
          channel_offset = mdh(idx).ushChannelId;
        end;

        mdh(idx).ushChannelId = mdh(idx).ushChannelId - channel_offset;

      end;

    end;  % IF "IS__VB13_PRE_RELEASE_VERSION"

    % finally, after MDH, read in the data
    adc = single(fread(fp, double(mdh(idx).ushSamplesInScan)*2, 'float32'));
    adc_real = adc(1:2:end);
    adc_imag = adc(2:2:end);
    adc_cplx = K_ICE_AMPL_SCALE_FACTOR * complex(adc_real,adc_imag);


    if ( DO__APPLY_FFT_SCALEFACTORS ),
      % scale whatever data comes in (QUESTION: should the noise data be scaled?)
      adc_cplx = adc_cplx * fft_scale( pos(03) );
    end;

    if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_REFLECT) ),
      if ( ~FLAG_data_reflect ),
        if ( DO__REVERSE_LINES ),
          disp('REFLECT detected, reversing lines.');
        else,
          disp('REFLECT detected.');
        end;
        FLAG_data_reflect = 1;
      end;
      if ( DO__REVERSE_LINES ),
        adc_cplx = flipud(adc_cplx);
      end;
    end;

    % for testing, abort after first repetition
    if ( mdh(idx).ushRepetition > 0 ),
      if ( ~FLAG_multiple_repetitions ),
        FLAG_multiple_repetitions = 1;
        disp('aborting after first repetition');
      end;

      % PHASCOR2D HACK: early versions of "epi_seg_3d" collect a dummy scan
      % for the 2D phascor after the phase correction reference lines and
      % increments the "repetition" loop counter for the 2D phascor lines
      if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), ...
                                  EvalInfoMask.MDH_PHASCOR) ),
        continue;
      elseif ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), ...
                                      EvalInfoMask.MDH_ONLINE) ),
        break;
      end;


    end;

    ACQEND = read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_ACQEND);
    if ( ACQEND ),


      % the last MDH and ADC line contains a 16-sample measurement of
      % unknown utility, so its stored in 'adc_cplx' without being
      % overwritten, but is not returned.
      break;

    end;  % IF "ACQEND"


    %------------------------------------------------------------------------%
    % store measurement lines

    pos = [
        1, ...                           % pos(01), columns
        1 + mdh(idx).ushLine, ...        % pos(02)
        1 + mdh(idx).ushChannelId, ...   % pos(03)
        1 + mdh(idx).ushSet, ...         % pos(04)
        1 + mdh(idx).ushEcho, ...        % pos(05)
        1 + mdh(idx).ushPhase, ...       % pos(06)
        1 + mdh(idx).ushRepetition, ...  % pos(07)
        1 + mdh(idx).ushSeg, ...         % pos(08)
        1 + mdh(idx).ushPartition, ...   % pos(09)
        1 + mdh(idx).ushSlice, ...       % pos(10)
        1 + mdh(idx).ushIda, ...         % pos(11)
        1 + mdh(idx).ushIdb, ...         % pos(12)
        1 + mdh(idx).ushIdc, ...         % pos(13)
        1 + mdh(idx).ushIdd, ...         % pos(14)
        1 + mdh(idx).ushIde, ...         % pos(15)
        1 + mdh(idx).ushAcquisition ...  % pos(16)  % note: acquisition is same as average
          ];

    if ( DO__CANONICAL_REORDER_COIL_CHANNELS ),
      pos(03)  = coil_order(pos(03));
    end;


    % correct for idiosyncratic, noncontiguous indexing schemes employed by
    % Siemens [found in both VB11A and VB13A so far]

    % Acquisition: the last of reference scan batch [three are collected for
    % each slice (2D) or partition (3D)] is assigned an acquisition index
    % of 1
    % [TODO: fix this hack]
    if ( (dimensions.NAveMeas == 1) && (mdh(idx).ushAcquisition > 0) ),
      pos(16) = 1;
    end;

    % Segment: all lines acquired in the negative readout direction (labeled
    % with the MDH_REFLECT flag) are assigned to their own segment via the
    % loopcounter to be compatible with the "OnlineTSE" phase correction
    % functor. as a result, a multishot segmented acquisition with
    % 'NSegMeas' legitimate segments will have 2*NSegMeas segments in the
    % loop counter---the even segments will be normal lines and the odd
    % segments will be reflected lines (0-based indexing). we must then
    % collapse the even and odd lines from each shot into a single segment.
    % [TODO: fix this hack...not all sequences use this scheme]

    % first, check if this line is a PHASCOR line (in case we're at the first line)
    if ( ~FLAG_phascor1d && read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PHASCOR) ),
      FLAG_phascor1d = 1;
    end;

    % ASSUME that if any phascor lines have been collected that all lines
    % (e.g., data and "patrefscan" lines) follow the same convention with
    % the segment loop counter---but this may not be true!
    if ( FLAG_phascor1d && FLAG_data_reflect ),
      pos(08) = floor(pos(08)/2);
    end;


    %%% begin logic to determine category of each ADC line

    %%% CATEGORY #1: noise adjustment scan lines
    if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_NOISEADJSCAN) ),

      if ( ~FLAG_noiseadjscan ),
        disp('NOISEADJSCAN detected.');
        FLAG_noiseadjscan = 1;
        noise_line = 0;
      end;

      % increment noise adjust scan line counter after all channels' data is in
      if ( mdh(idx).ushChannelId == 0 ),
        noise_line = noise_line + 1;
      end;

      % TODO: establish line counter for noise scans
      noiseadjscan(:, noise_line, pos(03)) = adc_cplx;
      continue;

    %%% CATEGORY #2: iPAT ACS lines reference scan
    elseif ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PATREFSCAN) ),

      if ( ~FLAG_patrefscan ),
        disp('PATREFSCAN detected.');
        FLAG_patrefscan = 1;
      end;

      %%% CATEGORY #2A: phase correction reference line for iPAT ACS lines reference scan
      if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PHASCOR) ),

        if ( ~FLAG_patrefscan_phascor ),
          disp('PATREFSCAN_PHASCOR detected.');
          FLAG_patrefscan_phascor = 1;
          ref_line = 0;
        end;

        % reset reference scan line counter at the beginning of each slice
        if ( (mdh(idx).ushChannelId == 0) && ...
             read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_FIRSTSCANINSLICE) ),
          ref_line = 0;
        end;


        % increment reference scan line counter after all channels' data is in
        if ( mdh(idx).ushChannelId == 0 ),
          ref_line = ref_line + 1;
        end;
	
	patrefscan_phascor(...
            :,   ref_line, ...     % only center line and center partition collected
            pos(03), pos(04), pos(05), pos(06), pos(07), pos(08),       1, ...
            pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;

      %%% CATEGORY #2B: lines used for both iPAT ACS and data
      elseif ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PATREFANDIMASCAN) ),
        % there are three PAT reconstruction modes: Integrated, Separate,
        % and Averaging. in 'Integrated' mode, the ACS lines are acquired
        % in-line with the data, so some serve both as data and ACS, whereas
        % in 'Separate' mode the ACS lines are acquired separately (e.g., as
        % in accelerated EPI). ('Averaging' mode is not supported
        % here....yet.)

        % for 'Integrated' mode, the ACS and data lines share the same
        % loopcounter numbering scheme, so the direct approach to storing
        % the k-space data would yield an inefficient storage where omitted
        % lines of k-space from the acceleration are stored as zeros! for
        % now, let's be lazy and waste some memory since this mode is
        % probably just for single-acquisition structural (e.g., 3D-MPRAGE)
        % images, so we don't have to be clever in storing the normal data
        % lines.

        % NOTE: in this case, the last skipped k-space lines will not be
        % present, so the iPAT data will contain R-1 fewer lines than the
        % corresponding non-accelerated data. we could just append those R-1
        % lines of zeros to the end, but that would require knowing the
        % value of R. hmmm... nah, skip it.

        if ( ~FLAG_patrefandimascan ),
          disp('PATREFANDIMASCAN detected.');
          FLAG_patrefandimascan = 1;
        end;

        patrefscan(...
            :, pos(02), ...
            pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
            pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;

        % hack, hack, hack
        if ( isempty(data) ),
          data = single(complex(zeros(size(adc_cplx))));
        end;

        % store ADC line in listed position within data volume (using 1-based indices)
        data(...
            :, pos(02), ...
            pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
            pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
            = adc_cplx;

      %%% CATEGORY #2C: just an ordinary ACS line
      else,
        patrefscan(...
            :, pos(02), ...
            pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
            pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;

      end;
      continue;

    %%% CATEGORY #3: phase stabilization scans
    elseif ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_REFPHASESTABSCAN) | ...
             read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PHASESTABSCAN) ),

      if ( ~FLAG_phasestabscan ),
        disp('PHASESTABSCAN detected.');
        FLAG_phasestabscan = 1;
      end;

      if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_REFPHASESTABSCAN) ),
	% count up dem echos!
        echo_num = echo_num + 1;

	% HACK: since the phase stability scan is always set as echo 0 in
        % the loopcounters, to avoid clobbering the true first echo for the
        % reference scans here we save the phase stability as echo N+1 where
        % N is the number of echos.
	
        refphasestabscan(...
            :, pos(02), ...
            pos(03), pos(04), echo_num, pos(06), pos(07), pos(08), pos(09), ...
            pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
            = adc_cplx;
      end;

      if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PHASESTABSCAN) ),
        phasestabscan(...
            :, pos(02), ...
            pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
            pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
            = adc_cplx;
      end;
      continue;

    %%% CATEGORY #4: data
    else,

      %%% CATEGORY #4A: phase correction reference line for data
      if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PHASCOR) ),

        if ( ~FLAG_data_phascor1d ),
          disp('PHASCOR (1D) detected.');
          FLAG_data_phascor1d = 1;

          % check if first phase correction reference line is reflected
          if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_REFLECT) ),
            % "ge_functionals" reflects the first reference line, so the
            % EVEN indexed reference lines will correspond to the ODD
            % indexed data lines. yuck!
            FLAG_data_phascor1d_orderswap = 1;
          end;

          ref_line = 0;
          ref_part = 1;
        end;

        % reset reference scan line counter at the beginning of each slice
        if ( (mdh(idx).ushChannelId == 0) && ...
             read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_FIRSTSCANINSLICE) ),
          ref_line = 0;
        end;


        % increment reference scan line counter after all channels' data is in
        if ( mdh(idx).ushChannelId == 0 ),
          ref_line = ref_line + 1;
        end;

        if ( FLAG_data_phascor1d_orderswap ),
          if ( ref_line == 2 ),
            data_phascor1d(...
                :, 2, ...
                pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), ref_part, ...
                pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = ...
                data_phascor1d(...
                    :, 1, ...
                    pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), ref_part, ...
                    pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16));

            data_phascor1d(...
                :, 1, ...
                pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), ref_part, ...
                pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;
            ref_line = ref_line + 1;
          end;
        end;

        data_phascor1d(...
            :, ref_line, ...     % only center line and center partition collected
            pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), ref_part, ...
            pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;


        % reference scan partition counter should NEVER be reset


        % PE line and PE partition incrementing depends on whether
        % acquisition is 2D or 3D. so can check whether line and partition
        % loop counter values on first line are equal to the center
        % line/partition of k-space.

        % increment reference scan partition counter at the end of each partition
        if ( ( (mdh(idx).ushChannelId + 1) == mdh(idx).ushUsedChannels ) && ...
             read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_LASTSCANINSLICE) ...
             && ( mdh(idx).ushPartition == mdh(idx).ushKSpaceCentrePartitionNo ) ...
             && ( mdh(idx).ushKSpaceCentrePartitionNo > 0 ) ...    % since always == 0 for 2D imaging
             ),
          ref_part = ref_part + 1;
        end;

        % QUESTION: for two partitions, is the center == 0 or == 1?


%%%/        if (  ( (mdh(idx).ushChannelId + 1) == mdh(idx).ushUsedChannels ) && ...
%%%/              ( (mdh(idx).ushAcquisition + 1) > 1 )  ),
%%%/
%%%/          % for some multi-shot scans (e.g., ge_functionals), get reference
%%%/          % lines with each shot
%%%/          ref_ = ref_ + 1;
%%%/        end;

      %%% CATEGORY #4B: field map lines for data
      % (PHASCOR2D HACK: currently the only signifier for 2D phascor data is
      % the absence of the MDH_ONLINE and MDH_PHASCOR flags)
      elseif ( ~read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_ONLINE) ),

        if ( ~FLAG_data_phascor2d ),
          disp('PHASCOR (2D) detected.');
          FLAG_data_phascor2d = 1;
          ref_line = 0;
          ref_part = 1;
        end;

        % TODO: when sequence code is fixed, this will be deleted!
        % reset reference scan line counter at the beginning of each slice
        if ( (mdh(idx).ushChannelId == 0) && ...
             read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_FIRSTSCANINSLICE) ),
          ref_line = 0;
        end;

        % TODO: when sequence code is fixed, this will be deleted!
        % increment reference scan line counter after all channels' data is in
        if ( mdh(idx).ushChannelId == 0 ),
          ref_line = ref_line + 1;
        end;


        % TODO: when sequence code is fixed, ref_line will be replaced by pos(02)!
        data_phascor2d(...
            :, ref_line, ...     % only center line and center partition collected
            pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), ref_part, ...
            pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;

        % reference scan partition counter should NEVER be reset

        % increment reference scan partition counter at the end of each partition
        % (ASSUME here that these reference lines ONLY appear in 3D sequences)
        if ( ( (mdh(idx).ushChannelId + 1) == mdh(idx).ushUsedChannels ) ...
             && read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), ...
                                       EvalInfoMask.MDH_LASTSCANINSLICE) ),
          ref_part = ref_part + 1;
        end;


      %%% CATEGORY #4B: data lines
      else,

        % hack, hack, hack
        if ( isempty(data) ),
          data = single(complex(zeros(size(adc_cplx))));
        end;

        % store ADC line in listed position within data volume (using 1-based indices)
        data(...
            :, pos(02), ...
            pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
            pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
            = adc_cplx;
      end;
    end; % IF (to determine category assignment)
  end;  % WHILE

  catch,

    caught = lasterror;
    status = dbstatus;

    if ( feof(fp) ),
      caught.message = sprintf('end of file encountered at scan number %d before ACQEND!  (debug: line %d)', ...
                               scan_num, caught.stack.line);
      error(caught);
    end;

    rethrow(caught);

  end;

  fclose(fp);


  t1 = clock;
  runtime_seconds = etime(t1,t0);

  TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                 fix(runtime_seconds/60/60), ...
                 rem(fix(runtime_seconds/60), 60), ...
                 rem(fix(runtime_seconds), 60));

  dstr = sprintf('total read time = %s', TIME);
  disp(sprintf('<t> [%s]: %s', mfilename, dstr));


  if ( nnz(data) == 0 ),
    if ( nargout > 0 ),
      error('no data found!!!');
    else,
      disp('no data found!!!');
      return;
    end;
  end;


  disp(sprintf('data memory allocated:    %10.2f MB', getfield(whos('data'), 'bytes')/2^20));

  density = nnz(data) / numel(data);
  % sparsity = 1 - density;
  redundancy = 1 / density;

  if ( redundancy > 1 && density < 0.99 ),
    disp(sprintf('data sparsity detected---potential redundancy factor = [[ %5.1fX ]]', ...
                 redundancy));
  end;


  % if auxiliary data is not present, set arrays to empty arrays in case
  % user requested arrays as output arguments

  if ( ~FLAG_data_phascor1d ),
    data_phascor1d            = [];
  end;

  if ( ~FLAG_data_phascor2d ),
    data_phascor2d           = [];
  end;

  if ( ~FLAG_noiseadjscan ),
    noiseadjscan            = [];
  end;

  if ( ~FLAG_patrefscan ),
    patrefscan              = [];
  end;

  if ( ~FLAG_patrefscan_phascor ),
    patrefscan_phascor      = [];
  end;


  if ( nargout > 0 ),
    varargout{1} = data;
    varargout{2} = data_phascor1d;
    varargout{3} = data_phascor2d;
    varargout{4} = noiseadjscan;
    varargout{5} = patrefscan;
    varargout{6} = patrefscan_phascor;
  end;


  return;


%**************************************************************************%
function bit = read_meas__extract_bit(mdh_eval_info, FLAG_EvalInfoMask)

  bit = bitget(mdh_eval_info, FLAG_EvalInfoMask+1);

  return;


%**************************************************************************%
function flag = read_meas__extract_flag(mdh_eval_info, flag_str)

%  persistent EvalInfoMask
  EvalInfoMask = read_meas__evaluation_info_mask_definition;

  field_str = sprintf('MDH_%s', flag_str);
  if ( ~isfield(EvalInfoMask, field_str) ),
    error('flag %s is not a valid EvalInfoMask flag', flag_str);
  end;

  flag = bitget(mdh_eval_info, EvalInfoMask.(field_str)+1);


  return;


%**************************************************************************%
function mdh_struct = read_meas__mdh_struct(varargin)
% snarfed from <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h>

  mdh_struct = struct(...
           'ulFlagsAndDMALength', [], ...
                      'lMeasUID', [], ...
                 'ulScanCounter', [], ...
                   'ulTimeStamp', [], ...
                'ulPMUTimeStamp', [], ...
               'aulEvalInfoMask', [], ...
                 'sEvalInfoMask', [], ...
              'ushSamplesInScan', [], ...
               'ushUsedChannels', [], ...
                       'ushLine', [], ...
                'ushAcquisition', [], ...   % note: acquisition is same as average
                      'ushSlice', [], ...
                  'ushPartition', [], ...
                       'ushEcho', [], ...
                      'ushPhase', [], ...
                 'ushRepetition', [], ...
                        'ushSet', [], ...
                        'ushSeg', [], ...
                        'ushIda', [], ...
                        'ushIdb', [], ...
                        'ushIdc', [], ...
                        'ushIdd', [], ...
                        'ushIde', [], ...
                        'ushPre', [], ...
                       'ushPost', [], ...
         'ushKSpaceCentreColumn', [], ...
                      'ushDummy', [], ...
             'fReadOutOffcentre', [], ...
             'ulTimeSinceLastRF', [], ...
         'ushKSpaceCentreLineNo', [], ...
    'ushKSpaceCentrePartitionNo', [], ...
            'aushIceProgramPara', [], ...
                  'aushFreePara', [], ...
                         'flSag', [], ...
                         'flCor', [], ...
                         'flTra', [], ...
                 'aflQuaternion', [], ...
                  'ushChannelId', [], ...
                 'ushPTABPosNeg', []);

  return;


%**************************************************************************%
function EvalInfoMask = read_meas__evaluation_info_mask_definition(varargin)
% snarfed from <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/MdhProxy.h>
% //##ModelId=3AFAAF7801CF

  EvalInfoMask.MDH_ACQEND            = 0;
  EvalInfoMask.MDH_RTFEEDBACK        = 1;
  EvalInfoMask.MDH_HPFEEDBACK        = 2;
  EvalInfoMask.MDH_ONLINE            = 3;
  EvalInfoMask.MDH_OFFLINE           = 4;
  EvalInfoMask.MDH_SYNCDATA          = 5;   % readout contains synchronous data
  EvalInfoMask.six                   = 6;
  EvalInfoMask.seven                 = 7;
  EvalInfoMask.MDH_LASTSCANINCONCAT  = 8;   % Flag for last scan in concatenation
  EvalInfoMask.nine                  = 9;

  EvalInfoMask.MDH_RAWDATACORRECTION = 10;  % Correct the rawdata with the rawdata correction factor
  EvalInfoMask.MDH_LASTSCANINMEAS    = 11;  % Flag for last scan in measurement
  EvalInfoMask.MDH_SCANSCALEFACTOR   = 12;  % Flag for scan specific additional scale factor
  EvalInfoMask.MDH_2NDHADAMARPULSE   = 13;  % 2nd RF excitation of HADAMAR
  EvalInfoMask.MDH_REFPHASESTABSCAN  = 14;  % reference phase stabilization scan
  EvalInfoMask.MDH_PHASESTABSCAN     = 15;  % phase stabilization scan
  EvalInfoMask.MDH_D3FFT             = 16;  % execute 3D FFT
  EvalInfoMask.MDH_SIGNREV           = 17;  % sign reversal
  EvalInfoMask.MDH_PHASEFFT          = 18;  % execute phase fft
  EvalInfoMask.MDH_SWAPPED           = 19;  % swapped phase/readout direction
  EvalInfoMask.MDH_POSTSHAREDLINE    = 20;  % shared line
  EvalInfoMask.MDH_PHASCOR           = 21;  % phase correction data
  EvalInfoMask.MDH_PATREFSCAN        = 22;  % additional scan for PAT reference line/partition
  EvalInfoMask.MDH_PATREFANDIMASCAN  = 23;  % additional scan for PAT reference line/partition that is also used as image scan
  EvalInfoMask.MDH_REFLECT           = 24;  % reflect line
  EvalInfoMask.MDH_NOISEADJSCAN      = 25;  % noise adjust scan --> Not used in NUM4
  EvalInfoMask.MDH_SHARENOW          = 26;  % all lines are acquired from the actual and previous e.g. phases
  EvalInfoMask.MDH_LASTMEASUREDLINE  = 27;  % indicates that the current line is the last measured line of all succeeding e.g. phases
  EvalInfoMask.MDH_FIRSTSCANINSLICE  = 28;  % indicates first scan in slice (needed for time stamps)
  EvalInfoMask.MDH_LASTSCANINSLICE   = 29;  % indicates last scan in slice (needed for time stamps)
  EvalInfoMask.MDH_TREFFECTIVEBEGIN  = 30;  % indicates the begin time stamp for TReff (triggered measurement)
  EvalInfoMask.MDH_TREFFECTIVEEND    = 31;  % indicates the end time stamp for TReff (triggered measurement)

  EvalInfoMask.thirty_two            = 32;
  EvalInfoMask.thirty_three          = 33;
  EvalInfoMask.thirty_four           = 34;
  EvalInfoMask.thirty_five           = 35;
  EvalInfoMask.thirty_six            = 36;
  EvalInfoMask.thirty_seven          = 37;
  EvalInfoMask.thirty_eight          = 38;
  EvalInfoMask.thirty_nine           = 39;

  EvalInfoMask.MDH_FIRST_SCAN_IN_BLADE       = 40;  % Marks the first line of a blade
  EvalInfoMask.MDH_LAST_SCAN_IN_BLADE        = 41;  % Marks the last line of a blade
  EvalInfoMask.MDH_LAST_BLADE_IN_TR          = 42;  % Set for all lines of the last BLADE in each TR interval


  EvalInfoMask.MDH_RETRO_LASTPHASE           = 45;  % Marks the last phase in a heartbeat
  EvalInfoMask.MDH_RETRO_ENDOFMEAS           = 46;  % Marks an ADC at the end of the measurement
  EvalInfoMask.MDH_RETRO_REPEATTHISHEARTBEAT = 47;  % Repeat the current heartbeat when this bit is found
  EvalInfoMask.MDH_RETRO_REPEATPREVHEARTBEAT = 48;  % Repeat the previous heartbeat when this bit is found
  EvalInfoMask.MDH_RETRO_ABORTSCANNOW        = 49;  % Just abort everything
  EvalInfoMask.MDH_RETRO_LASTHEARTBEAT       = 50;  % This adc is from the last heartbeat (a dummy)
  EvalInfoMask.MDH_RETRO_DUMMYSCAN           = 51;  % This adc is just a dummy scan, throw it away
  EvalInfoMask.MDH_RETRO_ARRDETDISABLED      = 52;  % Disable all arrhythmia detection when this bit is found


  return;


% VB13A ICE manual

%% There are 16 ICE Dimensions listed in the table below:
%% Dimension   Description
%% COL         Column related to pixel property ICE_PP_FIX
%% LIN         Line related to pixel property ICE_PP_FIX
%% CHA         Element related to pixel property ICE_PP_FIX
%% SET         Set related to pixel property ICE_PP_FIX
%% ECO         Contrast (echo) related to pixel property ICE_PP_VAR
%% PHS         Phase related to pixel property ICE_PP_FIX
%% REP         Measurement Repeat (Repetition) related to pixel property ICE_PP_VAR
%% SEG         Segment related to pixel property ICE_PP_FIX
%% PAR         Partition related to pixel property ICE_PP_VAR
%% SLC         Slice related to pixel property ICE_PP_VAR
%% IDA         1st free dispose of an ICE Dimension related to pixel property ICE_PP_VAR
%% IDB         2nd free dispose of an ICE Dimension related to pixel property ICE_PP_VAR
%% IDC         3rd free dispose of an ICE Dimension related to pixel property ICE_PP_VAR
%% IDD         4th free dispose of an ICE Dimension related to pixel property ICE_PP_FIX
%% IDE         5th free dispose of an ICE Dimension related to pixel property ICE_PP_FIX
%% AVE         Average related to pixel property ICE_PP_VAR
%%
%% Table 14: 16 ICE Dimensions in the order of generation in the IceProF.
%% There are 11 IceDimension with an intended use (e.g. SLC for slices)
%% and further 5 IceDimensions which can be used freely.


  %************************************************************************%
  %%% $Source: /space/repo/1/dev/dev/matlab/read_meas_dat.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
