function [d t MDHTime MPCUTime ecg2] = ReadSiemensPhysio(fname)
% [d t MDHTime MPCUTime ecg2] = ReadSiemensPhysio(fname)

% Data are 12 bit, so 4095 maximum. Values outside that range are
% some sort of trigger information or header.

% 6002 = begining of data (ECG Only)
% 5003 = end of data
% 5000 = Trigger on (rising edge, relevant)

% First 4 values are dont-care
% For Pulse and Resp, data immediately follow, end with 5003
% For ECG, there are 4 more values of unknown purpose, 
%          followed by 5002,
%          followed by header
%          followed by 6002
%          followed by data (VBx - 2 channels interleaved)
%          followed by 5003
%   Channel 1: values between 0 and 4095 (2048 is 0)
%   Channel aVF: values between 8192 and 12287 (10240 is 0)
% For Pulse, values of 5000 have something to do with triggering
% 

% Sample frequencies: ECG 400 Hz, Pulse 50 Hz, RESP 50 Hz.

% Is there a way to auto-detect the type of file? Eg, first 4 values.
% Pulse: 1 2 40 280
% Resp:  1 2 20 2 
% ECG:   1 1 2 40 

%fname = 'PhysioLog_20080805T161619.ecg';
%fname = 'PhysioLog_20080805T161619.puls';

d = [];
t = [];
ecg2 = [];
MDHTime  = [];
MPCUTime = [];

if(nargin ~= 1)
  fprintf('[d t MDHTime MPCUTime ecg2] = ReadSiemensPhysio(fname)\n')
  return;
end

fp = fopen(fname,'r');
if(fp == -1)
  fprintf('ERROR: opening %s\n',fname);
  return;
end

IsPulseFile = 0;
IsRespFile = 0;
IsECGFile = 0;

% Read 1st 4 values
h = fscanf(fp,'%d',4);
if(h(3) == 40) 
  IsPulseFile = 1;
  Tsamp = 1/50;
  fprintf('Detected as Pulse File\n');
elseif(h(3) == 20) 
  IsRespFile = 1;
  Tsamp = 1/50;
  fprintf('Detected as Resp File\n');
elseif(h(3) ==  2) 
  IsECGFile = 1;
  Tsamp = 1/400;
  fprintf('Detected as ECG File\n');
else
  fprintf('h = %d %d %d, Unrecognized\n',h);
  return;
end

if(IsECGFile)
  % Get to start of data
  s = '0';
  while(1 & ~isempty(s))
    s = fscanf(fp,'%s',1);
    if(strcmp(s,'6002')) break; end
  end
end


d = zeros(100000,1);
if(IsECGFile) ecg2 = zeros(100000,1); end
nth = 1;
while(1)
  s = fscanf(fp,'%s',1);
  if(strcmp(s,'5003')) break; end
  d(nth) = sscanf(s,'%d',1);
  if(IsECGFile)
    s = fscanf(fp,'%s',1);
    ecg2(nth) = sscanf(s,'%d',1);
  end
  nth = nth + 1;
end
d = d(1:nth-1);
if(IsECGFile) ecg2 = ecg2(1:nth-1); end
t = Tsamp*[0:length(d)-1];

s = ' '; 
while(1)
  s = fscanf(fp,'%s',1);
  if(isempty(s)) break; end
  if(strcmp(s,'LogStartMDHTime:'))
    LogStartMDHTime = fscanf(fp,'%d',1);
  elseif(strcmp(s,'LogStopMDHTime:'))
    LogStopMDHTime = fscanf(fp,'%d',1);
  elseif(strcmp(s,'LogStartMPCUTime:'))
    LogStartMPCUTime = fscanf(fp,'%d',1);
  elseif(strcmp(s,'LogStopMPCUTime:'))
    LogStopMPCUTime = fscanf(fp,'%d',1);
  end
end

MDHTime =  LogStopMDHTime - LogStartMDHTime;
MPCUTime = LogStopMPCUTime - LogStartMPCUTime;

fclose(fp);


