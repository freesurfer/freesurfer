from samseg.command_arguments import parse_sbtiv_args
from samseg.run_utilities import run_sbtiv_from_cmdargs

def sbtiv_main():
    run_sbtiv_from_cmdargs(parse_sbtiv_args())