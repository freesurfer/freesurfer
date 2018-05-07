import sys
import re
import platform


# check if file descriptor is a terminal device
def istty():
  return hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()


# a simple logging class that can be used to send all standard output to a file
class Logger(object):
  def __init__(self, filename):
    self.terminal = sys.stdout
    self.filename = filename
    open(self.filename, 'w').close()  # overwrite previous
    self.tty = istty()

  def write(self, message):
    if self.tty:
      self.terminal.write(message)
    with open(self.filename, 'a') as f:
      f.write(self.removeColor(message))

  def flush(self):
    self.terminal.flush()

  def removeColor(self, colorful):
    colorless = re.sub('\\033\[\d\dm', '', colorful)
    colorless = re.sub('\\033\[\dm', '', colorless)
    return colorless


# a class for colorful terminal output. Colors should be used by
# calling the global term object defined below
class TerminalColors(dict):
  def __init__(self):
    self.black = self.setCode('\033[30m')
    self.red = self.setCode('\033[31m')
    self.green = self.setCode('\033[32m')
    self.yellow = self.setCode('\033[33m')
    self.blue = self.setCode('\033[34m')
    self.magenta = self.setCode('\033[35m')
    self.cyan = self.setCode('\033[36m')
    self.white = self.setCode('\033[37m')
    self.underline = self.setCode('\033[4m')
    self.bold = self.setCode('\033[1m')
    self.end = self.setCode('\033[0m')

  def setCode(self, code):
    if istty() and platform.system() in ("Darwin", "Linux"): return code
    else: return ''


# create a 'global' TerminalColors
term = TerminalColors()


# print a warning message
def warning(message):
  print(term.yellow + '[warning] '+ term.end + message)


# print an error message
def error(message):
  print(term.red + '[error] ' + term.end + message)


# print an error message and exit
def errorExit(message, ret=1):
  error(message)
  sys.exit(ret)
