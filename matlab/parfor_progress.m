function fname = parfor_progress(fname,N)
% parfor_progress
%
% Displays the progress percent of the estimation.
%
% Adapted from Jeremy Scheff
% jdscheff@gmail.com - http://www.jeremyscheff.com/

if strcmp(fname,'init')
    i = 1;
    while exist(['parfor_progress' num2str(i) '.txt'], 'file') == 2
        i = i + 1;
    end;
    fname = ['parfor_progress' num2str(i) '.txt'];
    f = fopen(fname, 'w');
    if f < 0
        error('Do you have write permissions for %s?', pwd);
    end
    fprintf(f, '%d\n', N); % Save N at the top of file fname
    fclose(f);
else
    if N > 0
        f = fopen(fname, 'a');
        fprintf(f, [num2str(N) '\n']);
        fclose(f);
        f = fopen(fname, 'r');
        progress = fscanf(f, '%d');
        fclose(f);
        percent = 100*sum(progress(2:end))/progress(1);
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        disp([ 'Aproximate percentage of fitted locations: '  perc]);
    else
        delete(fname);
        disp('Aproximate percentage of fitted locations: 100%');    
    end;
end;
    

