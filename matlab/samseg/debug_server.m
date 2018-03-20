% USAGE: When at a breakpoint in MATLAB simply run this script from
% the command window. It will then look for {varname}.request files in the
% MATLAB_DUMP_DIR and write out {varname}.mat files with the variable data.
% The debugger server on Python side can then consume this data.
% To quit the debug server simply press ctrl-c once. This will return execution
% to the script you were running so you can continue to step or continue execution
% in the main script.

matlabDumpDir = getenv('MATLAB_DUMP_DIR')
while true
    listing = dir([matlabDumpDir '/*.request']);
    if length(listing) > 0
        for i = 1:length(listing)
            filename = listing(i).name;
            varname = filename(1:length(filename)-length('request')-1);
            if( eval(['exist(''' varname ''')']) == 1)
                disp(['Wrote ' varname])
                save([matlabDumpDir '/' varname '.mat'], varname)
            else
                disp(['Wrote ' varname ' as undefined'])
                eval([varname ' = ''undefined''']);
                save([matlabDumpDir '/' varname '.mat'], varname);
                eval(['clear ' varname]);
            end
            delete([matlabDumpDir '/' filename]);
        end
    end
    pause(0.1)
end
