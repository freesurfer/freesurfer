import sys
import re
import platform
import argparse
import textwrap


def istty():
    """Checks if file descriptor is a terminal device."""
    return hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()


class _TerminalColors(dict):
    """A class for colorful terminal output."""


    def __init__(self):
        self.black = self._setCode('\033[30m')
        self.red = self._setCode('\033[31m')
        self.green = self._setCode('\033[32m')
        self.yellow = self._setCode('\033[33m')
        self.blue = self._setCode('\033[34m')
        self.magenta = self._setCode('\033[35m')
        self.cyan = self._setCode('\033[36m')
        self.white = self._setCode('\033[37m')
        self.underline = self._setCode('\033[4m')
        self.bold = self._setCode('\033[1m')
        self.dim = self._setCode('\033[2m')
        self.end = self._setCode('\033[0m')


    def _setCode(self, code):
        if istty() and platform.system() in ("Darwin", "Linux"): return code
        else: return ''


# create a global, preinitialized TerminalColors
term = _TerminalColors()


def warning(message):
    """Prints a warning message."""
    print(term.yellow + 'warning: '+ term.end + message)


def error(message):
    """Prints an error message."""
    print(term.red + 'error: ' + term.end + message)


def errorExit(message, ret=1):
    """Prints an error message and exits with specified return code."""
    error(message)
    sys.exit(ret)


# global verbosity levels
_verbosity = 0


def setVerbosity(level):
    """Sets the global verbosity level."""
    global _verbosity
    _verbosity = level


def verbose(message, level=1):
    """Prints message dependent on the global verbosity level."""
    if _verbosity >= level:
        print(message)


# ----------------- Argument Parsing and Help Formatting -----------------


class _HelpFormatter(argparse.HelpFormatter):
    """HelpFormatter subclass used to format the help output to the freesurfer standard."""


    # Set default indent to 4 spaces
    def __init__(self, indent_increment=4, *args, **kwargs):
        super(_HelpFormatter, self).__init__(indent_increment=indent_increment, *args, **kwargs)


    # Format each argument
    def _format_action(self, action):
        parts = [self._current_indent * ' ' + self._format_action_invocation(action)]
        if action.help:
            self._indent()
            help_text = self._expand_help(action)
            for line in self._split_lines(help_text, self._width - self._current_indent):
                parts.append(self._current_indent * ' ' + line)
            self._dedent()
        return '\n'.join(parts) + '\n\n'


    # Format the argument invocation syntax
    def _format_action_invocation(self, action):
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return term.bold + metavar + term.end
        else:
            msg = ', '.join([term.bold + opt + term.end for opt in action.option_strings])
            if action.nargs != 0:
                default = self._get_default_metavar_for_optional(action)
                msg += '  ' +  term.blue + self._format_args(action, default) + term.end
            return msg


    # Format the usage string
    def _format_usage(self, usage, actions, groups, prefix):
        if prefix is None:
            prefix = 'usage: '

        # If usage is specified, use that
        if usage is not None:
            usage = usage % dict(prog=self._prog)

        # If no optionals or positionals are available, usage is just prog
        elif usage is None and not actions:
            usage = '%(prog)s' % dict(prog=self._prog)

        # If optionals and positionals are available, calculate usage
        elif usage is None:
            prog = '%(prog)s' % dict(prog=self._prog)

            # Split flags from positionals
            positionals = [action for action in actions if not action.option_strings]
            flags = [action for action in actions if action.option_strings]

            # Build full usage string
            format = self._format_actions_usage
            action_usage = format(positionals + flags, groups)
            usage = ' '.join([s for s in [prog, action_usage] if s])

            # Wrap the usage parts if it's too long
            text_width = self._width - self._current_indent
            if action_usage and len(prefix) + len(usage) > text_width:
                usage = prog + ' '
                first_len = len(prefix) + len(usage) + self._current_indent
                lines = self._split_lines(action_usage, self._width - first_len)
                usage += lines[0]
                for line in lines[1:]:
                    usage += '\n' + first_len * ' ' + line

        return self._current_indent * ' ' + prefix + usage + '\n\n'


    # Format the usage string arguments
    def _format_actions_usage(self, actions, groups):
        parts = []
        for i, action in enumerate(actions):
            if not action.required:
                continue
            if not action.option_strings:
                parts.append(self._format_args(action, self._get_default_metavar_for_positional(action)))
            else:
                parts.append(action.option_strings[0])
                if action.nargs != 0:
                    default = self._get_default_metavar_for_optional(action)
                    parts.append(self._format_args(action, default))
        parts.append('[options]')
        return ' '.join(parts)


    # Use type for flagged argument metavar
    def _get_default_metavar_for_optional(self, action):
        if action.type:
            return action.type.__name__
        else:
            return action.dest


    # Don't remove newlines when adding text
    def _fill_text(self, text, width, indent):
        return '\n'.join([textwrap.fill(line, width, initial_indent=indent, subsequent_indent=indent) for line in text.strip().splitlines()])


    # Don't remove newlines when adding help text to arguments
    def _split_lines(self, text, width):
        lines = []
        for line in text.strip().splitlines():
            lines.extend(textwrap.fill(line, width).splitlines())
        return lines


    class _Section(argparse.HelpFormatter._Section):


        # Reformat each section's header
        def format_help(self):
            # Format the indented section
            if self.parent is not None:
                self.formatter._indent()
            join = self.formatter._join_parts
            item_help = join([func(*args) for func, args in self.items])
            if self.parent is not None:
                self.formatter._dedent()

            # Return nothing if the section was empty
            if not item_help:
                return ''

            # Add the heading if the section was non-empty
            if self.heading is not argparse.SUPPRESS and self.heading is not None:
                indent = self.formatter._current_indent * ' '
                heading = indent + term.bold + self.heading + term.end + '\n'
            else:
                heading = ''

            # Join the section-initial newline, the heading and the help
            return join(['\n', heading, item_help, '\n'])


class ArgParser(argparse.ArgumentParser):
    """ArgumentParser subclass to format the help interface to the freesurfer standard."""


    # Use the custom freesurfer help formatter by default
    def __init__(self, add_help=True, allow_abbrev=False, formatter_class=_HelpFormatter, *args, **kwargs):
        super(ArgParser, self).__init__(allow_abbrev=allow_abbrev, formatter_class=formatter_class, add_help=False, *args, **kwargs)
        # Let's organize default arguments by required/optional status instead of positional/flagged
        self._action_groups.clear()
        self._required = self.add_argument_group('REQUIRED ARGUMENTS')
        self._optionals = self.add_argument_group('OPTIONAL ARGUMENTS')
        # _positionals is used as a default group in some places, so we'll remap it to _required
        self._positionals = self._required
        # By default, argparse adds the help flag at the beginning of the argument list. To get the help
        # flag at the end, we'll add it manually in parse_args().
        self.add_help = add_help


    # Wrap the default parse_args
    def parse_args(self, args=None, namespace=None):
        # Print usage and exit if no inputs provided.
        if not args and len(sys.argv[1:]) == 0:
            self.print_usage()
            self.exit(2)
        # Add help flag (at the end instead of the beginning)
        if self.add_help:
            self.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
        # Parse
        return super(ArgParser, self).parse_args(args, namespace)


    # Make sure default args are added to the appropriate custom groups
    def _add_action(self, action):
        if action.required:
            self._required._add_action(action)
        else:
            self._optionals._add_action(action)
        return action


    # Format the help text layout
    def format_help(self):
        formatter = self._get_formatter()

        # Usage
        formatter.start_section('USAGE')
        formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups, prefix='')
        formatter.end_section()

        # Description
        if self.description:
            formatter.start_section('DESCRIPTION')
            formatter.add_text(self.description)
            formatter.end_section()

        # Positionals, optionals and user-defined groups
        for action_group in self._action_groups:
            formatter.start_section(action_group.title)
            formatter.add_text(action_group.description)
            formatter.add_arguments(action_group._group_actions)
            formatter.end_section()

        # Epilog
        formatter.add_text(self.epilog)

        # Determine help from format above
        return formatter.format_help()
