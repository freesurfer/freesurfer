#!/usr/bin/env python3

import argparse
import logging
import os
import re
import sys

# logging.basicConfig(level="INFO")
_logger = logging.getLogger(__name__)


def build_argument_parser():
    desc = """Process C/C++ files for '#if 0' and '#if 1' blocks"""

    epilog = """This tool is for removing dead code from C/C++ files.
    It searches for '#if 0' and '#if 1' preprocessor directives (along
    with their '#else' clauses), and removes them along with any code
    which they disable. All other '#if', '#ifdef' and '#ifundef'
    directives are left untouched.
    """

    parser = argparse.ArgumentParser(description=desc, epilog=epilog)
    parser.add_argument("--input-file-name", help="Path to the input file",
                        required=True)
    parser.add_argument("--output-file-name", help="Path to the output file",
                        required=True)

    return parser


# -----------------------------------------------------

patterns = {
    "if0": re.compile(r"^\s*#if\s+0\s+"),
    "if1": re.compile(r"^\s*#if\s+1\s+"),
    "ifother": re.compile(r"^\s*#if\s+[^ 0-1]\w*"),
    "ifdef": re.compile(r"^\s*#ifdef\s+[^ 0-1]\w*"),
    "ifndef": re.compile(r"^\s*#ifndef\s+[^ 0-1]\w*"),
    "else": re.compile(r"^\s*#else\s+"),
    "endif": re.compile(r"^\s*#endif"),
    }

_START_KEYS = [ "if0", "if1", "ifother", "ifdef", "ifndef" ]
_START_PRESERVE_KEYS = [ "ifother", "ifdef", "ifndef" ]
_ELSE_KEY = "else"
_END_KEY = "endif"
_PLAIN_KEY = "plain"

def discover_line_type(line):
    result = _PLAIN_KEY
    match_names = []
    _logger.debug("Processing: %s", line)
    for name, pattern in patterns.items():
        _logger.debug("Considering %s", name)
        if pattern.match(line):
            _logger.debug("Matched %s", name)
            match_names.append(name)
    assert len(match_names) < 2
    if len(match_names)==1:
        result = match_names[0]
    return result

_PRESERVE_KEY = "preserve"
_ACTIVE_KEY = "active"

def check_active(if_stack):
    if len(if_stack)==0:
        return True

    checks = [x[_ACTIVE_KEY] for x in if_stack]
    return all(checks)
        

def create_if_stack_item(line_type):
    """Creates a dictionary for the if_stack

    The dictionary has two keys. One indicates whether the
    current #if block is active or not (i.e. whether the
    contents should go into the output). For our purposes,
    the only "inactive" start is "#if 0".

    The second key indicates whether or not the preprocessor
    directives should be preserved. We are only eliminating
    "#if 0" and "#if 1" blocks. Any "#else" and "#endif"
    statements belonging to other "#if" or "#ifdef" statements
    must be preseved.
    """
    assert line_type in _START_KEYS

    result = dict()
    result[_ACTIVE_KEY] = line_type != "if0"
    result[_PRESERVE_KEY] = line_type in _START_PRESERVE_KEYS

    return result


def process_lines(lines):
    _logger.info("Processing lines")
    result_lines = []

    # A stack of the various #if statements currently relevant
    if_stack = []
    for l in lines:
        pop_if_stack = False
        on_preproc_element = True
        line_type = discover_line_type(l)

        # See if the line is a preprocessor directive we care about
        if line_type in _START_KEYS:
            # We're starting a new block so push the stack
            if_stack.append(create_if_stack_item(line_type))
        elif line_type == _ELSE_KEY:
            # An '#else' will apply to the last #if
            last_element = if_stack[-1]
            if not last_element[_PRESERVE_KEY]:
                # We only flip whether we are active or not
                # for "#if 0" and "#if 1" blocks
                last_element[_ACTIVE_KEY] = not last_element[_ACTIVE_KEY]
            if_stack[-1] = last_element
        elif line_type == _END_KEY:
            # Flag that we've ended a block
            # Have to use this because we might need to output the "#endif"
            pop_if_stack = True
        elif line_type == _PLAIN_KEY:
            on_preproc_element = False
        else:
            raise ValueError("Unrecognised line_type: {0}".format(line_type))

        append_current_line = False
        if check_active(if_stack):
            if on_preproc_element:
                if if_stack[-1][_PRESERVE_KEY]:
                    append_current_line = True
            else:
                append_current_line = True
        if append_current_line:
            result_lines.append(l)
                    
        # Remove the last element in the stack if we just exited a block
        if pop_if_stack:
            if_stack.pop()

    # If the file is well-formed and we processed it correctly then
    # We should not be inside any preprocessor directives
    assert len(if_stack) == 0
        
    _logger.info("Line processing complete");
    return result_lines


def process_file(input_file_name, output_file_name):
    _logger.info("Reading file {}".format(input_file_name))
    text_lines = []
    with open(input_file_name, 'r') as input_file:
        text_lines = input_file.readlines()

    result_lines = process_lines(text_lines)
        
    _logger.info("Writing file {}".format(output_file_name))
    with open(output_file_name, 'w') as output_file:
        output_file.writelines(result_lines)

# -------------------------------------------------------
        
def main(argv):
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    process_file(args.input_file_name, args.output_file_name)


if __name__ == "__main__":
    main(sys.argv[1:])
