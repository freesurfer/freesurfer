function indfreq = indfreqfft(freq,Ntp,TR)
% indfreq = indfreqfft(freq,Ntp,TR)

if(nargin ~= 3) 
  msg = 'USAGE: indfreq = indfreqfft(freq,Ntp,TR)';
  qoe(msg); error(msg);
end

df = 1/(Ntp*TR); % frequency increment

fr = freq/df + 1; % ratio of frequency to increment

if( rem(fr,1) ~= 0 )
  % indfreq = [floor(fr) ceil(fr)];
  indfreq = round(fr);
else
  indfreq = fr;
end

return;
