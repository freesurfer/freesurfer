function S = readMetronic(CSVfile)

[num,txt,raw] = xlsread(CSVfile);

for n = 1:size(raw,1)
    if isnumeric(raw{n,1})
        nHdr = n-2;
        break
    end
end
S.hdr = raw(1:nHdr,:);
S.fields = raw(nHdr+1,:);
S.rows = raw(nHdr+2:end,:);
