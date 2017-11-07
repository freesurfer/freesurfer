def create_cmdargs(varargin):
    # mfileversion = '$Id$';
    #
    # % This is a structure to handle reading in of command-line args.  When
    # % adding a new arg, create a new field with the default value, then
    # % add a case in parse_args(). If the arg needs to be checked, add a
    # % case in check_params. Maybe add something to dump_args to print out
    # % the arg value. For clarity, do not reference this structure
    # % outside of this mfile.
    # cmdargs.involfiles = '';
    # cmdargs.missingStructures = '';
    # cmdargs.outdir = '';
    # cmdargs.regmatfile = '';
    # cmdargs.nthreads = 1;
    # cmdargs.debug = 0;
    # cmdargs.exvivo = 0;
    # %% Print useage if there are no arguments %%
    # if(nargin == 0)
    #   print_usage(cmdargs)
    #   return;
    # end
    # %% Parse the command-line arguments %%
    # cmdargs = parse_args(cmdargs,varargin);
    # if(isempty(cmdargs)) return; end
    # cmdargs = check_params(cmdargs);
    # if(isempty(cmdargs)) return; end
    # dump_args(cmdargs);
    #
    # fprintf('%s\n',mfileversion);
    # fprintf('Matlab version %s\n',version);
    pass


def parse_args(cmdargs, varargin):
    # %--------------------------------------------------%
    # % ----------- Parse Input Arguments ---------------%
    # function cmdargs = parse_args(cmdargs,varargin)
    #
    # inputargs = varargin{1};
    # ninputargs = length(inputargs);
    #
    # narg = 1;
    # while(narg <= ninputargs)
    #
    #   flag = deblank(inputargs{narg});
    #   narg = narg + 1;
    #   if(cmdargs.debug) fprintf(1,'Argument: %s\n',flag); end
    #   if(~isstr(flag))
    #     flag
    #     fprintf(1,'ERROR: All Arguments must be a string\n');
    #     error;
    #   end
    #
    #   switch(flag)
    #
    #    case '--i',
    #     arg1check(flag,narg,ninputargs);
    #     cmdargs.involfiles = strvcat(cmdargs.involfiles,inputargs{narg});
    #     narg = narg + 1;
    #
    #    case '--threads',
    #     arg1check(flag,narg,ninputargs);
    #     cmdargs.nthreads = sscanf(inputargs{narg},'%d');
    #     narg = narg + 1;
    #
    #    case '--o',
    #     arg1check(flag,narg,ninputargs);
    #     cmdargs.outdir = inputargs{narg};
    #     narg = narg + 1;
    #
    #    case '--regmat',
    #     arg1check(flag,narg,ninputargs);
    #     cmdargs.regmatfile = inputargs{narg};
    #     narg = narg + 1;
    #
    #    case '--missing',
    #     arg1check(flag,narg,ninputargs);
    #     cmdargs.missingStructures = strvcat(cmdargs.missingStructures,inputargs{narg});
    #     narg = narg + 1;
    #
    #    case '--debug',
    #     cmdargs.debug = 1;
    #
    #    case '--exvivo',
    #     cmdargs.exvivo = 1;
    #
    #    otherwise
    #     fprintf(2,'ERROR: Flag %s unrecognized\n',flag);
    #     cmdargs = [];
    #     return;
    #
    #   end % --- switch(flag) ----- %
    #
    # end % while(narg <= ninputargs)
    #
    # return;
    # %--------------------------------------------------%
    return cmdargs


def check_params(cmdargs):
    # %--------------------------------------------------%
    # %% Check argument consistency, etc %%%
    # function cmdargs = check_params(cmdargs)
    #   if(isempty(cmdargs.involfiles))
    #     fprintf('ERROR: must specify at least one input with --i\n');
    #     error;
    #   end
    #   if(isempty(cmdargs.outdir))
    #     fprintf('ERROR: must specify an output dir with --o\n');
    #     error;
    #   end
    # return;
    return cmdargs


def arg1check(flag, nflag, nmax):
    # %--------------------------------------------------%
    # %% Check that there is at least one more argument %%
    # function arg1check(flag,nflag,nmax)
    #   if(nflag>nmax)
    #     fprintf(1,'ERROR: Flag %s needs one argument',flag);
    #     error;
    #   end
    # return;
    pass


def arg2check(flag, nflag, nmax):
    # %--------------------------------------------------%
    # %% Check that there are at least two more arguments %%
    # function arg2check(flag,nflag,nmax)
    #   if(nflag > nmax-1 )
    #     fprintf(1,'ERROR: Flag %s needs two arguments',flag);
    #     error;
    #   end
    # return;
    pass


def print_usage(cmdargs):
    # %------------- Print Usage ---------------------%
    # function print_usage(cmdargs)
    #   fprintf('USAGE:\n');
    #   fprintf('run_samseg\n');
    #   fprintf(' --o <outdir>          : output folder\n');
    #   fprintf(' --i <input>           : input volume (add --i input for more)\n');
    #   fprintf(' --threads <nthreads>  : number of threads\n');
    #   fprintf(' --regmat <regmat>     : regmat from previous run\n');
    #   fprintf(' --missing <struct>    : specify a missing structure\n');
    #   fprintf(' --exvivo              : run samseg exvivo');
    # return
    pass


def dump_args(cmd_args):
    # %------------- Dump Args ---------------------%
    # function dump_args(cmdargs)
    #   fprintf(' outdir %s\n',cmdargs.outdir);
    #   for n = 1:size(cmdargs.involfiles)
    #     fprintf(' %d %s\n',n,cmdargs.involfiles(n,:));
    #   end
    #   fprintf('nthreads %d\n',cmdargs.nthreads);
    #   if(~isempty(cmdargs.regmatfile))
    #     fprintf('regmatfile %s\n',cmdargs.regmatfile)
    #   end
    #
    # return
    pass


