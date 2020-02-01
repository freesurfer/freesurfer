import sys
import argparse
import textwrap

from . import term


class HelpFormatter(argparse.HelpFormatter):

    # set default indent to 4 spaces
    def __init__(self, indent_increment=4, *args, **kwargs):
        super(HelpFormatter, self).__init__(indent_increment=indent_increment, *args, **kwargs)

    # format each argument
    def _format_action(self, action):
        parts = [self._current_indent * ' ' + self._format_action_invocation(action)]
        if action.help:
            self._indent()
            help_text = self._expand_help(action)
            for line in self._split_lines(help_text, self._width - self._current_indent):
                parts.append(self._current_indent * ' ' + line)
            self._dedent()
        return '\n'.join(parts) + '\n\n'

    # format the argument invocation syntax
    def _format_action_invocation(self, action):
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return term.bold(metavar)
        else:
            msg = ', '.join([term.bold(opt) for opt in action.option_strings])
            if action.nargs != 0:
                default = self._get_default_metavar_for_optional(action)
                msg += '  ' +  term.blue(self._format_args(action, default).upper())
            return msg

    # format the usage string
    def _format_usage(self, usage, actions, groups, prefix):
        if prefix is None:
            prefix = 'usage: '

        # if usage is specified, use that
        if usage is not None:
            usage = usage % dict(prog=self._prog)

        # if no optionals or positionals are available, usage is just prog
        elif usage is None and not actions:
            usage = '%(prog)s' % dict(prog=self._prog)

        # if optionals and positionals are available, calculate usage
        elif usage is None:
            prog = '%(prog)s' % dict(prog=self._prog)

            # split flags from positionals
            positionals = [action for action in actions if not action.option_strings]
            flags = [action for action in actions if action.option_strings]

            # build full usage string
            format = self._format_actions_usage
            action_usage = format(positionals + flags, groups)
            usage = ' '.join([s for s in [prog, action_usage] if s])

            # wrap the usage parts if it's too long
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

    # use type for flagged argument metavar
    def _get_default_metavar_for_optional(self, action):
        if action.type:
            return action.type.__name__
        else:
            return action.dest

    # don't remove newlines when adding text
    def _fill_text(self, text, width, indent):
        return '\n'.join([textwrap.fill(line, width, initial_indent=indent, subsequent_indent=indent) for line in text.strip().splitlines()])

    # don't remove newlines when adding help text to arguments
    def _split_lines(self, text, width):
        lines = []
        for line in text.strip().splitlines():
            lines.extend(textwrap.fill(line, width).splitlines())
        return lines


    class _Section(argparse.HelpFormatter._Section):

        # reformat each section's header
        def format_help(self):
            # format the indented section
            if self.parent is not None:
                self.formatter._indent()
            join = self.formatter._join_parts
            item_help = join([func(*args) for func, args in self.items])
            if self.parent is not None:
                self.formatter._dedent()

            # return nothing if the section was empty
            if not item_help:
                return ''

            # add the heading if the section was non-empty
            if self.heading is not argparse.SUPPRESS and self.heading is not None:
                indent = self.formatter._current_indent * ' '
                heading = indent + term.bold(self.heading) + '\n'
            else:
                heading = ''

            # join the section-initial newline, the heading and the help
            return join(['\n', heading, item_help, '\n'])


class ExtendAction(argparse.Action):
    '''
    Custom action to extend arguments into a list so that both of these
    are allowed:

    command --input 1 2 3
    command --input 1 --input 2 --input 3
    '''
    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest) or []
        items.extend(values)
        setattr(namespace, self.dest, items)


class ArgumentParser(argparse.ArgumentParser):
    '''ArgumentParser subclass to format the help interface to the freesurfer standard.'''

    # use the custom freesurfer help formatter by default
    def __init__(self, add_help=True, allow_abbrev=False, formatter_class=HelpFormatter, *args, **kwargs):
        super(ArgumentParser, self).__init__(allow_abbrev=allow_abbrev, formatter_class=formatter_class, add_help=False, *args, **kwargs)
        # let's organize default arguments by required/optional status instead of positional/flagged
        self._action_groups.clear()
        self._required = self.add_argument_group('REQUIRED ARGUMENTS')
        self._optionals = self.add_argument_group('OPTIONAL ARGUMENTS')
        # _positionals is used as a default group in some places, so we'll remap it to _required
        self._positionals = self._required
        # By default, argparse adds the help flag at the beginning of the argument list. To get the help
        # flag at the end, we'll add it manually in parse_args()
        self.add_help = add_help
        # register our custom extend action
        self.register('action', 'extend', ExtendAction)

    # wrap the default parse_args
    def parse_args(self, args=None, namespace=None):
        # print usage and exit if no inputs are provided
        if not args and not sys.argv[1:] and self._required._group_actions:
            self.print_usage()
            self.exit(2)
        # add help flag (at the end instead of the beginning)
        if self.add_help:
            self.add_argument('-h', '--help', action='help', help='Show this help message and exit.')
        return super(ArgumentParser, self).parse_args(args, namespace)

    # make sure default args are added to the appropriate custom groups
    def _add_action(self, action):
        if action.required:
            self._required._add_action(action)
        else:
            self._optionals._add_action(action)
        return action

    # format the help text layout
    def format_help(self):
        formatter = self._get_formatter()

        # usage
        formatter.start_section('USAGE')
        formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups, prefix='')
        formatter.end_section()

        # description
        if self.description:
            formatter.start_section('DESCRIPTION')
            formatter.add_text(self.description)
            formatter.end_section()

        # positionals, optionals and user-defined groups
        for action_group in self._action_groups:
            formatter.start_section(action_group.title)
            formatter.add_text(action_group.description)
            formatter.add_arguments(action_group._group_actions)
            formatter.end_section()

        # epilog
        formatter.add_text(self.epilog)

        # determine help from format above
        return formatter.format_help()
