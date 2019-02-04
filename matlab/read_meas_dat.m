function varargout = read_meas_dat(filename, user_options)
%READ_MEAS_DAT  read in Siemens format raw data, VB13A/VB15A-style "meas.dat"
%
% data = read_meas_dat(filename, <options>)
%
% [data, phascor1d, phascor2d, noise] = read_meas_dat(filename)
%
% [data, phascor1d, phascor2d, noise, patrefscan, patrefscan_phascor] = read_meas_dat(filename)
%
% [data, phascor1d, phascor2d, noise, patrefscan, patrefscan_phascor, phasestabscan, refphasestabscan] = read_meas_dat(filename)
%
% "options" is a structure which allows changing of some of the defaults:
%
%    options.ReverseLines                 -- set to 0 or 1  (default is 1)
%    options.PhascorCollapseSegments      -- set to 0 or 1  (default is 0)
%    options.CanonicalReorderCoilChannels -- set to 0 or 1  (default is 1)
%    options.ApplyFFTScaleFactors         -- set to 0 or 1  (default is 0)
%    options.ReadMultipleRepetitions      -- set to 0 or 1  (default is 1)
%    options.ReturnStruct                 -- set to 0 or 1  (default is 0)
%    options.MatrixDouble                 -- set to 0 or 1  (default is 0)
%
% if a field does not exist, then the default is used. "options" is optional :).
%
%
% the raw k-space data is returned in a 16-dimensional array, following
% the convention used in ICE. the mapping between the array dimensions
% and the loopcounters used in ICE is as follows:
%
%    #01:  ColMeas
%    #02:  LinMeas
%    #03:  ChaMeas
%    #04:  SetMeas
%    #05:  EcoMeas
%
%    #06:  PhsMeas
%    #07:  RepMeas
%    #08:  SegMeas
%    #09:  ParMeas
%    #10:  SlcMeas
%
%    #11:  IdaMeas
%    #12:  IdbMeas
%    #13:  IdcMeas
%    #14:  IddMeas
%    #15:  IdeMeas
%
%    #16:  AveMeas


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

% 2007/jul/08: added fix for non-contiguous coil channels found in 7T data

% 2007/jul/09: improved support for "meas.out" files lacking headers

% 2007/aug/07: added ability to return all data as fields of a struct

% 2007/aug/14: fixed phase stabilization support for multiecho acquisitions

% 2007/oct/16: fixed order swapping for even and odd phascor lines in EPI

% 2007/oct/19: added ability to return data even if file incomplete


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 10/04/2006
% $Id: read_meas_dat.m,v 1.12 2008/07/16 06:45:24 jonnyreb Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.12 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %------------------------------------------------------------------------%
  % basic error checking

  matlab_version = version;
  if ( str2num(matlab_version(1:3)) < 7.2 )
    disp(sprintf('"%s" only supported for MATLAB 7.2 and higher', mfilename));
    return;
  end;

  if ( ~exist(filename, 'file') ),
    error('file [%s] does not exist', filename);
  end;


  %------------------------------------------------------------------------%
  % the pre-release version installed on the Bay4 128-channel host in 2005
  % differs slightly from the full-release version, and files must be read
  % in differently depending on the version. yuck!
  IS__VB13_PRE_RELEASE_VERSION = 0;
  IS__VB13A_RELEASE_VERSION = 1;

  DO__MDH_SAVE = 0;

  DO__RECOVER_FROM_INCOMPLETE = 1;

  
  
  DO__FLIP_REFLECTED_LINES = 1;
  DO__PHASCOR_COLLAPSE_SEGMENTS = 0;
  DO__CANONICAL_REORDER_COIL_CHANNELS = 0;
  DO__APPLY_FFT_SCALEFACTORS = 0;
  DO__READ_MULTIPLE_REPETITIONS = 1;

  % by default, return individual arrays (for backward compatibility)
  DO__RETURN_STRUCT = 0;

  % by default, return all data as full arrays of type SINGLE
  DO__MATRIX_DOUBLE = 0;
  DO__MATRIX_SPARSE = 0;


  % if the "options" struct is provided by caller, then override default
  % flags

  if ( exist('user_options', 'var') ),

    if ( isfield(user_options, 'ReverseLines') ),
      DO__FLIP_REFLECTED_LINES = user_options.ReverseLines;
      disp(sprintf(' :FLIP_REFLECTED_LINES = %d', DO__FLIP_REFLECTED_LINES));
    end;

    if ( isfield(user_options, 'PhascorCollapseSegments') ),
      DO__PHASCOR_COLLAPSE_SEGMENTS = user_options.PhascorCollapseSegments;
      disp(sprintf(' :PHASCOR_COLLAPSE_SEGMENTS = %d', DO__PHASCOR_COLLAPSE_SEGMENTS));
    end;

    if ( isfield(user_options, 'CanonicalReorderCoilChannels') ),
      DO__CANONICAL_REORDER_COIL_CHANNELS = user_options.CanonicalReorderCoilChannels;
      disp(sprintf(' :CANONICAL_REORDER_COIL_CHANNELS = %d', DO__CANONICAL_REORDER_COIL_CHANNELS));
    end;

    if ( isfield(user_options, 'ApplyFFTScaleFactors') ),
      DO__APPLY_FFT_SCALEFACTORS = user_options.ApplyFFTScaleFactors;
      disp(sprintf(' :APPLY_FFT_SCALEFACTORS = %d', DO__APPLY_FFT_SCALEFACTORS));
    end;

    if ( isfield(user_options, 'ReadMultipleRepetitions') ),
      DO__READ_MULTIPLE_REPETITIONS = user_options.ReadMultipleRepetitions;
      disp(sprintf(' :READ_MULTIPLE_REPETITIONS = %d', DO__READ_MULTIPLE_REPETITIONS));
    end;

    if ( isfield(user_options, 'ReturnStruct') ),
      DO__RETURN_STRUCT = user_options.ReturnStruct;
      disp(sprintf(' :RETURN_STRUCT = %d', DO__RETURN_STRUCT));
    end;

    if ( isfield(user_options, 'MatrixDouble') ),
      DO__MATRIX_DOUBLE = user_options.MatrixDouble;
      disp(sprintf(' :MATRIX_DOUBLE = %d', DO__MATRIX_DOUBLE));
    end;

  end;

  options = struct('ReverseLines',                  DO__FLIP_REFLECTED_LINES, ...
                   'PhascorCollapseSegments',       DO__PHASCOR_COLLAPSE_SEGMENTS, ...
                   'CanonicalReorderCoilChannels',  DO__CANONICAL_REORDER_COIL_CHANNELS, ...
                   'ApplyFFTScaleFactors',          DO__APPLY_FFT_SCALEFACTORS, ...
                   'ReadMultipleRepetitions',       DO__READ_MULTIPLE_REPETITIONS, ...
                   'ReturnStruct',                  DO__RETURN_STRUCT, ...
                   'MatrixDouble',                  DO__MATRIX_DOUBLE);




  % constants defined in <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h>
  MDH_NUMBEROFEVALINFOMASK   = 2;
  MDH_NUMBEROFICEPROGRAMPARA = 4;

  MDH_FREEHDRPARA = 4;

  % from "MrServers/MrVista/include/Ice/IceDefs.h":
  ICE_RAWDATA_SCALE       = 131072.0;  % 64 ^ 3 / 2
  K_ICE_AMPL_SCALE_FACTOR = 80 * 20 * ICE_RAWDATA_SCALE / 65536;


  %------------------------------------------------------------------------%

  t0 = clock;

  [pathstr, filestr, extstr] = fileparts(filename);

  [fp, errstr] = fopen(filename, 'r', 'l');
  if ( fp == -1 ),
    error(errstr);
  end;


  % determine size (in bytes) of ascii header files stored in the 'meas.dat'
  % format (i.e., "Config_.evp", "Dicom_.evp", etc) to skip over them all.
  % [note: this works for VB11A also---the first integer is 32, which is
  % the number of bytes to skip at the beginning!]
  data_start = fread(fp, 1, 'uint32');

  % if file has no header (e.g., if its in "meas.out" format),
  % "data_start" should be 32.
  if ( data_start == 32 ),
    disp('no header detected! assuming "meas.out" format...');

    IS__VB13A_RELEASE_VERSION = 0;

    [pathstr, name, ext] = fileparts(filename);
    meas_asc = fullfile(pathstr, [name, '.asc']);

    if ( exist(meas_asc, 'file') ),
      disp('discovered corresponding "meas.asc" file! parsing...');

      [asc, errstr] = fopen(meas_asc, 'r', 'l');
      if ( asc == -1 ), error(errstr); end;

      % read header into one string for parsing
      header = fscanf(asc, '%c');

      fclose(asc);

    else,
      % set header to empty string
      header = '';
    end;

    % can't sort channels without header information  :(
    DO__CANONICAL_REORDER_COIL_CHANNELS = 0;

    % jump to beginning of binary data, let's get to work!
    fseek(fp, data_start, 'bof');

  else,

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

%%      % the number of segments is listed in two places in the header with the
%%      % field names "NSegMeas" and "NSeg", and for some reason only the "NSeg"
%%      % field gives the correct number of segments; SO for this field we break
%%      % with the convention [TODO: fix the bug in the martinos EPI sequence
%%      % that causes this hack]
%%      % UPDATE: it appears that "NSegMeas" corresponds to the cumulative
%%      % number of distinct segments appearing in the loop counters, whereas
%%      % "NSeg" is the true number of segments that is the same as the
%%      % number of shots, i.e., they should differ by a factor of 2 when the
%%      % "OnlineTSE" functor is in use.
%%      if ( strcmp(param, 'NSegMeas') ),
%%        param = 'NSeg';
%%      end;

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

    % VB11A workaround hack (for EPI only)
    if ( ~IS__VB13A_RELEASE_VERSION ),
      dimensions.NColMeas = 0;
      dimensions.NAveMeas = 1;
      dimensions.NSegMeas = 2;
    end;


    %------------------------------------------------------------------------%
    % extract FFT scalefactors (if utility file is in path)

    if ( DO__APPLY_FFT_SCALEFACTORS && ...
         exist('read_meas_dat__fft_scalefactors', 'file') ),
      [fft_scale, fft_scale_channel] = read_meas_dat__fft_scalefactors(header);
      if ( ~isempty( fft_scale ) ),
        disp('applying FFT scale factors.');
      else,
        DO__APPLY_FFT_SCALEFACTORS = 0;
        disp('FFT scale factors not found!!!');
      end;
    else,
      fft_scale = [];
      disp('ignoring FFT scale factors.');
    end;


    %------------------------------------------------------------------------%
    % compute canonical coil ordering from coil element strings---the channel
    % number assigned to each coil is determined by the order in which the
    % channels are selected in SYNGO, and NOT by the preamp channel numbers,
    % so to enforce a consistent ordering across scans that is independent of
    % the operator's coil selections we can sort by the fixed preamp channel
    % strings.

    if ( DO__CANONICAL_REORDER_COIL_CHANNELS && ...
         exist('read_meas_dat__reorder_coil_channels', 'file') ),
      disp('reordering coil channels to canonical preamp-based ordering.');

      [coil_index, coil_order] = read_meas_dat__reorder_coil_channels(header);

      % to map using index:
      %   dataval(coil_index);
      % to map using order:
      %   coil_order(dataind);

      if ( DO__APPLY_FFT_SCALEFACTORS ),
        if ( length(fft_scale) ~= length(coil_index) ),
          DO__APPLY_FFT_SCALEFACTORS = 0;
          disp('mismatch between number of channels and FFT scale factors!!!  -->  ignoring FFT scale factors...');
        else,
          % don't forget to reorder FFT scalefactors!
          fft_scale = fft_scale(coil_index);
        end;
      end;

    else,
      % clear flag if auxiliary function not in path
      DO__CANONICAL_REORDER_COIL_CHANNELS = 0;
    end;

    %------------------------------------------------------------------------%

    % jump past header to beginning of binary data
    fseek(fp, data_start, 'bof');

  end;  % end of header parsing when present, if..else..end

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

  mdh = read_meas__mdh_struct;

  %--------------------------------------------------------------------------%

  matrixtype = @single;
  matrixtypestr = 'single';

  % unfortunately, matlab does not support sparse arrays whose number of
  % dimensions is larger than 2.   :(
  if ( DO__MATRIX_SPARSE ),
    %%%matrixtype = @sparse;
  end;

  if ( DO__MATRIX_DOUBLE ),
    matrixtype = @double;
    matrixtypestr = 'double';
  end;


  % number of measurements along all 16 dimensions can be found in the
  % 'Config_.evp' file contained within "meas.dat" header; allocate memory
  % for measurement data, and grow data matrix with each line (this method
  % is, surprisingly, faster than pre-allocation!)
%  dimensions
%  tic;
%  fprintf(' allocating memory...');
%  data                    = complex(zeros(struct2array(dimensions), matrixtypestr));
%  fprintf('done!  ');
%  toc;
%  fprintf('\n');
  data                    = matrixtype([]);

  FLAG__data               = 0;
  FLAG__data_reflect       = 0;
  FLAG__data_swapped       = 0;

  FLAG__phascor1d          = 0;  % indicating presence of any type of phascor line

  data_phascor1d           = matrixtype([]);
  FLAG__data_phascor1d     = 0;
  data_phascor2d           = matrixtype([]);
%  scan_phascor2d           = matrixtype([]);
  FLAG__data_phascor2d     = 0;
  data_fieldmap            = matrixtype([]);
  FLAG__data_fieldmap      = 0;
  noiseadjscan             = matrixtype([]);
  FLAG__noiseadjscan       = 0;
  patrefscan               = matrixtype([]);
%  scan_patrefscan          = matrixtype([]);
  FLAG__patrefscan         = 0;
  patrefscan_phascor       = matrixtype([]);
  FLAG__patrefscan_phascor = 0;
  FLAG__patrefandimascan   = 0;
  phasestabscan            = matrixtype([]);
  refphasestabscan         = matrixtype([]);
  FLAG__phasestabscan      = 0;

  
  FLAG__data_phascor1d_orderswap = 0;

  FLAG__patrefscan_phascor_orderswap = 0;

  FLAG__syncdata = 0;

  % set to zero if trouble
  FLAG__status_OK = 1;


  % convert binary mask into strings
  EvalInfoMask = read_meas__evaluation_info_mask_definition;
  EvalInfoNames = fieldnames(EvalInfoMask);


  % channel index offset (needed only if IS__VB13_PRE_RELEASE_VERSION)
  channel_offset = NaN;

  % store all channel indices in case non-contiguous (e.g., 7T host)
  FLAG__stored_channel_indices = 0;
  channel_indices = [];

  FLAG__channel_remap = 0;
  channel_index_map = [];

  ACQEND = 0;

  try,
    %------------------------------------------------------------------------%

    ulDMALength = uint16(fread(fp, 1, 'uint16'));
    ulFlags1    = uint8(fread(fp, 1, 'uint8'));
    ulFlags2    = uint8(fread(fp, 1, 'uint8'));

    fseek(fp, fpos, 'bof');


    % loop through file, peeling off one MDH and ADC line at a time
    while ( ~ACQEND ),


      %----------------------------------------------------------------------%

      meas_pos = ftell(fp);
      meas_num = meas_num + 1;

      if ( DO__MDH_SAVE ),
        idx = meas_num;
      else,
        % overwrite mdh(1) with every measurement line
        idx = 1;
      end;

      ulDMALength = uint16(fread(fp, 1, 'uint16'));
%      fseek(fp, dimensions.NColMeas * 2 * 4 + 128 - 2, 'cof');
%      fseek(fp, meas_pos + ulDMALength*4, 'bof');
%      continue;

      % to skip this scan and start reading next MDH:
      %%% fseek(fp, ulDMALength-2, 'cof')
      %%% fseek(fp, meas_pos + ulDMALength, 'bof')

      % to slurp up MDH and data for this scan:
      %%% binary = fread(fp, double(ulDMALength-2), 'uchar=>uchar');

      ulFlags1    = uint8(fread(fp, 1, 'uint8'));
      ulFlags2    = uint8(fread(fp, 1, 'uint8'));

      % the old way, at beginning of each MDH:    mdh(idx).ulFlagsAndDMALength        = uint32(fread(fp, 1, 'uint32'));
      mdh(idx).ulFlagsAndDMALength        = uint32(double(ulFlags2) * 2^24) + uint32(double(ulFlags1) * 2^16) + uint32(ulDMALength);

      mdh(idx).lMeasUID                   = int32(fread(fp, 1, 'int32'));
      mdh(idx).ulScanCounter              = uint32(fread(fp, 1, 'uint32'));
      if ( ~isempty(mdh(idx).ulScanCounter) ), scan_num = mdh(idx).ulScanCounter; end;
      mdh(idx).ulTimeStamp                = uint32(round(fread(fp, 1, 'uint32')*2.5)); % milliseconds
      mdh(idx).ulPMUTimeStamp             = uint32(round(fread(fp, 1, 'uint32')*2.5)); % milliseconds

%      timestamp(mdh(idx).ulScanCounter).mdh = mdh(idx).ulTimeStamp;
%      timestamp(mdh(idx).ulScanCounter).pmu = mdh(idx).ulPMUTimeStamp;

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

      % for some auxiliary scans (e.g., "AdjCoilSens" for the prescan
      % normalize), the first scan contains gibberish and its mask contains
      % the 'SYNCDATA' bit. if this is found, just skip ahead to the next
      % MDH.
      if ( ~FLAG__syncdata && read_meas__extract_flag(mdh(idx).aulEvalInfoMask(1), 'SYNCDATA') ),
        fread(fp, double(ulDMALength) - 128, 'uint8');
        continue;
      end;
      FLAG__syncdata = 1;


      % finally, after MDH, read in the data
      adc = (fread(fp, double(mdh(idx).ushSamplesInScan)*2, 'float32=>single'));
      adc_real = adc(1:2:end);
      adc_imag = adc(2:2:end);
      adc_cplx = K_ICE_AMPL_SCALE_FACTOR * complex(adc_real,adc_imag);


      if ( DO__APPLY_FFT_SCALEFACTORS ),
        % scale whatever data comes in (QUESTION: should the noise data be scaled?)
        adc_cplx = adc_cplx * fft_scale( double(mdh(idx).ushChannelId+1) );
      end;

      if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_REFLECT) ),
        if ( ~FLAG__data_reflect ),
          if ( DO__FLIP_REFLECTED_LINES ),
            disp(' REFLECT detected, reversing lines.');
          else,
            disp(' REFLECT detected.');
          end;
          FLAG__data_reflect = 1;
        end;
        if ( DO__FLIP_REFLECTED_LINES ),
          adc_cplx = flipud(adc_cplx);
        end;
      end;

      if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_SWAPPED) ),
	if ( ~FLAG__data_swapped ),
	  % not sure what to do about swapped lines, but for now just report it to user.
	  disp(' SWAPPED detected, ignoring.');
          FLAG__data_swapped = 1;
	end;
      end;
	
      
      % for testing, abort after first repetition
      if ( ~DO__READ_MULTIPLE_REPETITIONS && mdh(idx).ushRepetition > 0 ),

%%%        % PHASCOR2D HACK: early versions of "epi_seg_3d" collect a dummy scan
%%%        % for the 2D phascor after the phase correction reference lines and
%%%        % increments the "repetition" loop counter for the 2D phascor lines
%%%        if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), ...
%%%                                    EvalInfoMask.MDH_PHASCOR) ),
%%%          continue;
%%%        elseif ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), ...
%%%                                        EvalInfoMask.MDH_ONLINE) ),
          disp('aborting after first repetition...');
          break;
%%%        end;
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


      % (all this is to cater to the 7T host's peculiarities. sigh...)
      if ( ~FLAG__stored_channel_indices ),
        if ( ~ismember(pos(03), channel_indices) ),
          % accumulate list of channel indices
          channel_indices(end+1) = pos(03);

          % record / update first and last channel numbers
          channel_1 = channel_indices(1) - 1;
          channel_N = channel_indices(end) - 1;
        else,
          % all channels have been stored
          FLAG__stored_channel_indices = 1;

          % weirdness found
          if ( (channel_indices(1) ~= 1) || ...
               (channel_indices(end) ~= length(channel_indices)) ),
            FLAG__channel_remap = 1;

            % the coils should at LEAST be in order even if they are not contiguous
            [channel_sort, channel_pos] = sort(channel_indices);
            channel_index_map = zeros(1,max(channel_indices));
            channel_index_map(channel_sort) = channel_pos;

            % after all the coil channels have been recorded, we need to go
            % back and resort the data that's been stored out-of-order.
            % since we don't know which data has been seen so far, test
            % everything.
            if ( ~isempty(data) ), data = data(:,:,channel_sort); end;

            if ( FLAG__data_phascor1d ),     data_phascor1d     = data_phascor1d(    :,:,channel_sort); end;
            if ( FLAG__data_phascor2d ),     data_phascor2d     = data_phascor2d(    :,:,channel_sort); end;
            if ( FLAG__data_fieldmap ),      data_fieldmap      = data_fieldmap(     :,:,channel_sort); end;
            if ( FLAG__noiseadjscan ),       noiseadjscan       = noiseadjscan(      :,:,channel_sort); end;
            if ( FLAG__patrefscan ),         patrefscan         = patrefscan(        :,:,channel_sort); end;
            if ( FLAG__patrefscan_phascor ), patrefscan_phascor = patrefscan_phascor(:,:,channel_sort); end;
            if ( FLAG__phasestabscan ),      refphasestabscan   = refphasestabscan(  :,:,channel_sort); end;

          end;
        end;
      end;

      % apply map
      if ( FLAG__channel_remap ),
        pos(03) = channel_index_map( pos(03) );
      end;

      if ( DO__CANONICAL_REORDER_COIL_CHANNELS ),
        if ( channel_1 ~= 0 ),
          disp('non-contiguous channels found---aborting canonical reordering.');
          DO__CANONICAL_REORDER_COIL_CHANNELS = 0;
        end;
        pos(03)  = coil_order(pos(03));
      end;


      % correct for idiosyncratic, noncontiguous indexing schemes employed by
      % Siemens [found in both VB11A and VB13A so far]

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
      if ( ~FLAG__phascor1d && read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PHASCOR) ),
        FLAG__phascor1d = 1;
      end;

      % ASSUME that if any phascor lines have been collected that all lines
      % (e.g., data and "patrefscan" lines) follow the same convention with
      % the segment loop counter---but this may not be true!
      if ( FLAG__phascor1d && FLAG__data_reflect && DO__PHASCOR_COLLAPSE_SEGMENTS ),
        pos(08) = pos(08)/2;  % if condition true, pos(08) should always be even, so quotient should be integer!
      end;

      % Acquisition: the last of reference navigators scan batch [three are
      % collected for each slice (2D) or partition (3D)] is assigned an
      % acquisition index of 1, for siemens's "EPIPhaseCorrPEFunctor"
      % [TODO: fix this hack]
      if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), ...
                                  EvalInfoMask.MDH_PHASCOR) ),
        pos(16) = ceil(pos(16)/2);
      end;


      %%% begin logic to determine category of each ADC line

      %%% CATEGORY #1: noise adjustment scan lines
      if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_NOISEADJSCAN) ),

        if ( ~FLAG__noiseadjscan ),
          disp(' NOISEADJSCAN detected.');
          FLAG__noiseadjscan = 1;
          noise_line = 0;
        end;

        % increment noise adjust scan line counter after all channels' data is in
        if ( mdh(idx).ushChannelId == channel_1 ),
          noise_line = noise_line + 1;
        end;

        % TODO: establish line counter for noise scans

        % remove fft scaling of noise
        if ( DO__APPLY_FFT_SCALEFACTORS ),
          adc_cplx = adc_cplx / fft_scale( double(mdh(idx).ushChannelId+1) );
        end;
        noiseadjscan(:, noise_line, pos(03)) = adc_cplx;

        continue;

        %%% CATEGORY #2: iPAT ACS lines reference scan
      elseif ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PATREFSCAN) ),

        if ( ~FLAG__patrefscan ),
          disp(' PATREFSCAN detected.');
          FLAG__patrefscan = 1;
        end;

        %%% CATEGORY #2A: phase correction navigator line for iPAT ACS lines reference scan
        if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PHASCOR) ),

          if ( ~FLAG__patrefscan_phascor ),
            disp(' PATREFSCAN_PHASCOR detected.');
            FLAG__patrefscan_phascor = 1;

            % check if first phase correction navigator line is reflected (see comments below)
            if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_REFLECT) ),
              FLAG__patrefscan_phascor_orderswap = 1;
            end;

	    nav_base_segment = pos(08);

            nav_line = 0;
	    
          end;

          % reset navigator scan line counter at the beginning of each slice
          if ( (mdh(idx).ushChannelId == channel_1) && ...
               read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_FIRSTSCANINSLICE) ),
            nav_line = 0;
          end;

          % reset navigator scan line counter
          if ( (mdh(idx).ushChannelId == channel_1) && ...
               pos(08) >= (nav_base_segment+2) )
            nav_line = 0;
	    nav_base_segment = pos(08);
          end;


          % increment navigator scan line counter after all channels' data is in
          if ( mdh(idx).ushChannelId == channel_1 ),
            nav_line = nav_line + 1;
          end;

          patrefscan_phascor(...
              :,   nav_line, ...     % only center line and center partition collected
              pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), 1, ...
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
          % lines of k-space are stored as zeros! for now, let's be lazy and
          % waste some memory since this mode is probably just for
          % single-acquisition structural (e.g., 3D-MPRAGE) images, so we
          % don't have to be clever in storing the normal data lines.

          % NOTE: in this case, the last skipped k-space lines will not be
          % present, so the iPAT data will contain R-1 fewer lines than the
          % corresponding non-accelerated data. we could just append those R-1
          % lines of zeros to the end, but that would require knowing the
          % value of R. hmmm... nah, skip it.

          % so if a line is both ACS and an image line, store it *twice*: once
          % as data and once as reference.

          % (for 1-dimensional acceleration, R could be estimated from the
          % data sparsity calculated at the end, but this would not work for
          % 2-dimensional imaging!)

          if ( ~FLAG__patrefandimascan ),
            disp(' PATREFANDIMASCAN detected.');
            FLAG__patrefandimascan = 1;
          end;

          patrefscan(...
              :, pos(02), ...
              pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
              pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;

          % hack, hack, hack
          if ( isempty(data) ),
            data = complex(zeros(size(adc_cplx), matrixtypestr));
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

%          scan_patrefscan(...
%              1, pos(02), ...
%              pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
%              pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = mdh(idx).ulScanCounter;

        end;
        continue;

        %%% CATEGORY #3: phase stabilization scans
      elseif ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PHASESTABSCAN) ),

        if ( ~FLAG__phasestabscan ),

          % need the time of the phase stabilization navigator echo, "TS",
          % for image reconstruction, but its NOT available from the header
          % for some weird reason. extract it from the MDH header. it should
          % be the same across segments and slices (which are the only
          % dimensions over which its repeated).
          TS = double(mdh(idx).ulTimeSinceLastRF) / 1e3;    % store as milliseconds

          disp(sprintf(' PHASESTABSCAN detected;  (( TS = %.3f ms )).', TS));
          FLAG__phasestabscan = 1;

        end;

        %%% (note to self) tested with the following examples:
        %%%  -  meas_MID864_gre_mgh_multiecho_FID16967.dat
        %%%  -  meas_MID865_gre_mgh_multiecho_body_FID16968.dat
        %%%  -  meas_MID1048_Flash_SWI_Best_8Ch_lowBW_FID6643.dat
        %%%  -  meas_MID1579_gre_2d_phasestab_FID4764.dat
        %%%  -  meas_MID1582_gre_3d_phasestab_FID4767.dat


        % NOTE: both the phase stabilization measurements and the phase
        % stabilization reference measurements are stored as echo 0
        % regardless of how many echoes are acquired for the imaging data.

        %%% VB13A ICE manual, V0.7, p. 274:
        %  "Phase stabilization is restricted to the sharing of phase
        %   stabilization scans between different contrasts (echoes), i.e.
        %   the imaging scans, which belong to the same phase stabilization
        %   echo, may only differ in the ECO ICE Dimension."

        if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_REFPHASESTABSCAN) ),

          % the phase stabilization works by collecting a full set
          % (potentially multiple echoes) of lines at the beginning of each
          % slice that are tagged as MDH_REFPHASESTABSCAN, followed by a
          % *single* measurement that is tagged as both MDH_REFPHASESTABSCAN
          % and MDH_PHASESTABSCAN. it is this last scan that serves as the
          % reference, and it will exhibit the same timing as the subsequent
          % phase stabilization scans that follow each group of image
          % echoes. therefore the group of MDH_REFPHASESTABSCAN lines
          % preceding each *single* reference echo are not used by the
          % correction and are discarded.

          refphasestabscan(...
              :, pos(02), ...
              pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
              pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
              = adc_cplx;
        else,

          phasestabscan(...
              :, pos(02), ...
              pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
              pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
              = adc_cplx;

        end;
        continue;

        %%% CATEGORY #4: data
      else,

        %%% CATEGORY #4A: phase correction navigator line for data
        if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_PHASCOR) ),

          if ( ~FLAG__data_phascor1d ),
            disp(' PHASCOR (1D) detected.');
            FLAG__data_phascor1d = 1;

            % check if first phase correction navigator line is reflected
            if ( read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_REFLECT) ),

              % for example, siemens's "ep2d_bold" reflects the first
              % navigator line, in which case the EVEN-indexed navigator
              % lines will correspond to the ODD-indexed data lines. yuck!
              FLAG__data_phascor1d_orderswap = 1;

              % TODO: store the reflected and non-reflected lines in
              % separate matrices, instead of adhering to the "odd--even"
              % organization, since siemens appears to like to collect more
              % of one than the other! (DONE!)
            end;

	    nav_base_segment = pos(08);
            
	    nav_line = 0;
	    
          end;

          % reset navigator scan line counter at the beginning of each slice
          if ( mdh(idx).ushChannelId == channel_1 && read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_FIRSTSCANINSLICE) ),
	    nav_line = 0;
	    nav_reset = nav_line;
	    nav_base_segment = pos(08);
          end;

          % reset navigator scan line counter
          if ( (mdh(idx).ushChannelId == channel_1) && ...
               pos(08) >= (nav_base_segment+2) )
            nav_line = nav_reset;
	    nav_base_segment = pos(08);
          end;


          % increment navigator scan line counter after all channels' data is in
          if ( mdh(idx).ushChannelId == channel_1 ),
            nav_line = nav_line + 1;
          end;
	  
	  % if 2D, only one partition. if 3D each navigator line is
	  % collected at center of k-space, but partition loopcounter set
	  % to center partition, so here we force it to 1 so that only
	  % nav_line counter is incremented.
	  
	  data_phascor1d(...
	      :, nav_line, ...     % only center line and center partition collected
	      pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), 1, ...
	      pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;
	  
	  
          % navigator scan partition counter should NEVER be reset


          % PE line and PE partition incrementing depends on whether
          % acquisition is 2D or 3D. so can check whether line and partition
          % loop counter values on first line are equal to the center
          % line/partition of k-space.

%          % increment navigator scan partition counter at the end of each partition
%          if ( FLAG__stored_channel_indices && ...
%               ( mdh(idx).ushChannelId == channel_N ) && ...
%               read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_LASTSCANINSLICE) ...
%               && ( mdh(idx).ushPartition == mdh(idx).ushKSpaceCentrePartitionNo ) ...
%               && ( mdh(idx).ushKSpaceCentrePartitionNo > 0 ) ...    % since always == 0 for 2D imaging
%               ),
%            nav_part = nav_part + 1;
%          end;

          % QUESTION: for two partitions, is the center == 0 or == 1?


          %%%/        if (  ( (mdh(idx).ushChannelId + 1) == mdh(idx).ushUsedChannels ) && ...
          %%%/              ( (mdh(idx).ushAcquisition + 1) > 1 )  ),
          %%%/
          %%%/          % for some multi-shot scans (e.g., ge_functionals), get nagivator
          %%%/          % lines with each shot
          %%%/          nav_ = nav_ + 1;
          %%%/        end;

          %%% CATEGORY #4B: field map lines for data
          % (PHASCOR2D HACK: currently the only signifier for 2D phascor data is
          % the absence of the MDH_ONLINE and MDH_PHASCOR flags)
        elseif ( ~read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), EvalInfoMask.MDH_ONLINE) ),

          if ( ~FLAG__data_phascor2d ),
            disp(' PHASCOR (2D) detected.');
            FLAG__data_phascor2d = 1;
            nav_part = 1;
          end;

          % navigator scan partition counter should NEVER be reset

          data_phascor2d(...
              :, pos(02), ...     
              pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), nav_part, ...  % only center partition collected
              pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;
	  
%	  scan_phascor2d(...
%              1, pos(02), ...     
%              pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), nav_part, ...  
%              pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = mdh(idx).ulScanCounter;
	  
	  
          % increment navigator scan partition counter at the end of each partition (signified with the "LASTSCANINSLICE" flag)
          % (ASSUME here that these navigator lines ONLY appear in 3D sequences)
          if ( FLAG__stored_channel_indices && ...
               ( mdh(idx).ushChannelId == channel_N ) ...
               && read_meas__extract_bit(mdh(idx).aulEvalInfoMask(1), ...
                                         EvalInfoMask.MDH_LASTSCANINSLICE) ),
            nav_part = nav_part + 1;
          end;
	  
	  
          %%% CATEGORY #4C: data lines (finally!)
        else,

          if ( ~FLAG__data ),
            FLAG__data = 1;
            fpos_data = ftell(fp);
          end;

          % hack, hack, hack
          if ( isempty(data) ),
            data = complex(zeros(size(adc_cplx), matrixtypestr));
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


    %------------------------------------------------------------------------%

  catch,

    caught = lasterror;
    status = dbstatus;

    disp(sprintf('<!> [%s]:  trouble at line %d...', mfilename, caught.stack(1).line));

    if ( feof(fp) ),

      caught.message = sprintf('EOF encountered before ACQEND! last scan number = %d ', ...
                               scan_num);

      if ( DO__RECOVER_FROM_INCOMPLETE ),

        fprintf(1, '\n');
        warning(caught.identifier, caught.message);
        fprintf(1, '\n');

        disp(sprintf('<!> [%s]:  ABORTING read, incomplete file!!!', mfilename));
        fprintf(1, '\n');

        FLAG__status_OK = 0;

      else,
        error(caught.identifier, caught.message);
      end;

    else,

      rethrow(caught);

    end;


  end;


  %------------------------------------------------------------------------%
  %%% post file I/O reorganization

  % if the phase correction lines and the data lines mismatch in terms
  % whether the odd or even lines are reflected, shift the phase
  % correction line index by one and insert a duplicate of the second
  % line at the beginning.


  if ( DO__PHASCOR_COLLAPSE_SEGMENTS && FLAG__patrefscan_phascor_orderswap ),

    disp('first "PATREFSCAN_PHASCOR" line reflected! compensating with odd--even order swap...');

    patrefscan_phascor(:, 2:end+1, :, :, :, :, :, :, :, :, :, :, :, :, :, :) = ...
        patrefscan_phascor(:, 1:end, :, :, :, :, :, :, :, :, :, :, :, :, :, :);

    patrefscan_phascor(:, 1, :, :, :, :, :, :, :, :, :, :, :, :, :, :) = ...
        patrefscan_phascor(:, 3, :, :, :, :, :, :, :, :, :, :, :, :, :, :);
  end;


  if ( DO__PHASCOR_COLLAPSE_SEGMENTS && FLAG__data_phascor1d_orderswap ),

    disp('first "PHASCOR" line reflected! compensating with odd--even order swap...');

    data_phascor1d(:, 2:end+1, :, :, :, :, :, :, :, :, :, :, :, :, :, :) = ...
        data_phascor1d(:, 1:end, :, :, :, :, :, :, :, :, :, :, :, :, :, :);

    data_phascor1d(:, 1, :, :, :, :, :, :, :, :, :, :, :, :, :, :) = ...
        data_phascor1d(:, 3, :, :, :, :, :, :, :, :, :, :, :,: , :, :);
  end;


  %------------------------------------------------------------------------%

  fclose(fp);


  t1 = clock;
  runtime_seconds = etime(t1,t0);

  TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                 fix(runtime_seconds/60/60), ...
                 rem(fix(runtime_seconds/60), 60), ...
                 rem(fix(runtime_seconds), 60));

  dstr = sprintf('total read time = %s', TIME);
  disp(sprintf('<t> [%s]: %s', mfilename, dstr));

%
%  if ( nnz(data) == 0 ),
%    if ( nargout > 0 ),
%      warning('no data found!!!');
%    else,
%      disp('no data found!!!');
%      return;
%    end;
%  end;
%
%
%  disp(sprintf('data memory allocated:    %10.2f MB', getfield(whos('data'), 'bytes')/2^20));
%
%  density = nnz(data) / numel(data);
%  % sparsity = 1 - density;
%  redundancy = 1 / density;
%
%  if ( redundancy > 1 && density < 0.99 ),
%    disp(sprintf('data sparsity detected---potential redundancy factor = [[ %5.1f X ]]', ...
%                 redundancy));
%  end;
%

  %------------------------------------------------------------------------%

  if ( DO__RETURN_STRUCT ),

    meas = struct;
    meas.file = strcat(filestr, extstr);
    meas.time = datestr(t1, 'yyyy-mmm-dd HH:MM:SS');

    meas.data = data;

    if ( FLAG__data_phascor1d ),
      meas.data_phascor1d = data_phascor1d;
    end;

    if ( FLAG__data_phascor2d ),
      meas.data_phascor2d = data_phascor2d;
%      meas.scan_phascor2d = scan_phascor2d;
    end;

    if ( FLAG__noiseadjscan ),
      meas.noiseadjscan = noiseadjscan;
    end;

    if ( FLAG__patrefscan ),
      meas.patrefscan = patrefscan;
%      meas.scan_patrefscan = scan_patrefscan;
    end;

    if ( FLAG__patrefscan_phascor ),
      meas.patrefscan_phascor = patrefscan_phascor;
    end;

    if ( FLAG__phasestabscan ),
      meas.phasestabscan = phasestabscan;
      meas.refphasestabscan = refphasestabscan;
    end;

    if ( DO__MDH_SAVE ),
      meas.mdh = mdh;
    end;

    if ( exist('read_meas_prot', 'file') && IS__VB13A_RELEASE_VERSION ),
      [meas.prot, meas.evp] = read_meas_prot(filename, header);
    end;

%    meas.timestamp = timestamp;

    meas.options = options;

    if ( FLAG__status_OK ),
      meas.STATUS = 'success';
    else,
      meas.STATUS = 'ABORTED';
    end;


    if ( nargout > 0 ),
      varargout{1} = meas;
    end;


  else,

    % if auxiliary data is not present, set arrays to empty arrays in case
    % user requested arrays as output arguments

    if ( ~FLAG__data_phascor1d ),
      data_phascor1d            = [];
    end;

    if ( ~FLAG__data_phascor2d ),
      data_phascor2d           = [];
    end;

    if ( ~FLAG__noiseadjscan ),
      noiseadjscan            = [];
    end;

    if ( ~FLAG__patrefscan ),
      patrefscan              = [];
    end;

    if ( ~FLAG__patrefscan_phascor ),
      patrefscan_phascor      = [];
    end;

    if ( ~FLAG__phasestabscan ),
      phasestabscan          = [];
      refphasestabscan       = [];
    end;


    if ( nargout > 0 ),
      varargout{1} = data;
      varargout{2} = data_phascor1d;
      varargout{3} = data_phascor2d;
      varargout{4} = noiseadjscan;
      varargout{5} = patrefscan;
      varargout{6} = patrefscan_phascor;
      varargout{7} = phasestabscan;
      varargout{8} = refphasestabscan;
    end;

  end;

  return;


%**************************************************************************%
function bit = read_meas__extract_bit(mdh_eval_info, FLAG__EvalInfoMask)

  bit = bitget(mdh_eval_info, FLAG__EvalInfoMask+1);

  return;


%**************************************************************************%
function flag = read_meas__extract_flag(mdh_eval_info, flag_str)

%slow!  persistent EvalInfoMask
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
