import sys
import re
import shutil
import platform


# terminal color codes
colors = {
    'black': '\033[30m',
    'red': '\033[31m',
    'green': '\033[32m',
    'yellow': '\033[33m',
    'blue': '\033[34m',
    'magenta': '\033[35m',
    'cyan': '\033[36m',
    'white': '\033[37m',
    'underline': '\033[4m',
    'dim': '\033[2m',
    'bold': '\033[1m',
    'end': '\033[0m',
}

# check if tty allows colors
istty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()
colorsEnabled = istty and platform.system() in ("Darwin", "Linux")
# clear all color codes if they're not enabled
if not colorsEnabled:
    colors = dict.fromkeys(colors, '')

# temporary function to dynamically create string-wrapping methods
def mkmethod(name, code):
    def wrapstring(string):
        return ''.join((str(code), string, colors['end']))
    return wrapstring

# create the color methods and remove the temporary variables
for color, code in colors.items():
    vars()[color] = mkmethod(color, code)
del color, code, mkmethod


def colorformat(string):
    '''Formats a string with the text formatting options available in `colors`.

    Args:
        string (str): String to format.
    Returns:
        The formatted string.
    '''
    return string.format(**colors)


def removecolor(string):
    '''Removes color codes from a string.

    Args:
        string (str): String to reformat.
    Returns:
        The string with color codes removed.
    '''
    colorless = re.sub('\\033\[\d\dm', '', string)
    colorless = re.sub('\\033\[\dm', '', colorless)
    return colorless


def width(fallback=60):
    '''Gets the current terminal width.

    Args:
        fallback (int): Default width if none is found. Defaults to 60.
    Returns:
        Terminal width.
    '''
    return shutil.get_terminal_size((fallback, 20)).columns
