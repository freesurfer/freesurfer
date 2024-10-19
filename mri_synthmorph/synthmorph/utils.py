import shutil
import textwrap


def resolve_abbrev(needle, strings, lower=False):
    """Return a full-length string matching a substring from the beginning.

    Parameters
    ----------
    needle : str
        Substring of one of several `strings`.
    strings : str or iterable of str
        Full-length strings, one of which should begin with `needle`.
    lower : bool, optional
        Convert needle to lowercase before matching.

    Returns
    -------
    str
        String in `strings` that begins with `needle` if there is no ambiguity.
        If there is not exactly one match, the function will return `needle`.

    """
    if isinstance(strings, str):
        strings = [strings]
    strings = tuple(strings)

    if lower:
        needle = needle.lower()

    matches = [f for f in strings if f.startswith(needle)]
    return matches[0] if len(matches) == 1 else needle


def rewrap_text(text, width=None, hard='\t\n', hard_indent=0, end=''):
    """Rewrap text such that lines fill the available horizontal space.

    Reformats individual paragraphs of a text body, considering subsequent
    lines with identical indentation as paragraphs. For unspecified width, the
    function will attempt to determine the extent of the terminal.

    Parameters
    ----------
    text : str
        Text to rewrap.
    width : int, optional
        Maximum line width. None means the width of the terminal as determined
        by `textwrap`, defaulting to 80 characters for background processes.
    hard : str, optional
        String interpreted as a hard break when terminating a line. Useful for
        inserting a line break without changing the indentation level. Must end
        with a line break and will be removed from the output text.
    hard_indent : int, optional
        Number of additional whitespace characters by which to indent the lines
        following a hard break. See `hard`.
    end : str, optional
        Append to the reformatted text.

    Returns
    -------
    out : str
        Reformatted text.

    """
    # Inputs.
    if width is None:
        width = shutil.get_terminal_size().columns
    lines = text.splitlines(keepends=True)

    # Merge lines to paragraphs.
    pad = []
    pad_hard = []
    par = []
    for i, line in enumerate(lines):
        ind = len(line) - len(line.lstrip())
        if i == 0 or ind != pad[-1] or lines[i - 1].endswith(hard):
            par.append('')
            pad.append(ind)
            pad_hard.append(ind)

        if line.endswith(hard):
            line = line.replace(hard, '\n')
            pad_hard[-1] += hard_indent
        par[-1] += line[ind:]

    # Reformat paragraphs.
    for i, _ in enumerate(par):
        par[i] = textwrap.fill(
            par[i], width,
            initial_indent=' ' * pad[i], subsequent_indent=' ' * pad_hard[i],
        )

    return '\n'.join(par) + end
