from samseg.command_arguments import parse_args
from samseg.run_samseg_ported import run_samseg_from_cmdargs

def samseg_main():
    run_samseg_from_cmdargs(parse_args())