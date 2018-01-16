import collections
import re
import csv
from strings import *

class Tokenizer:

    def __init__(self):
        self.Token = collections.namedtuple('Token', ['type', 'value', 'line', 'start', 'end'])

    def tokenize(self, text):
        token_specification = [
            ('HEADER',    r'=+\s*[a-zA-Z ]+\s*=+'),     # Header: start of a section.
            ('COLUMN',    r'[a-zA-Z ]+'),               # Column: data columns.
            ('DATA',      r'[0-9.]+'),                  # Data: data to process.
            ('NEWLINE',   r'\n'),                       # Newline character.
            ('MISC',      r'\w*'),                      # Blank lines and things.
            ('UNDEFINED', r'.'),                        # No token match, skip.
        ]
        tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)
        line_num = 1
        line_start = 0
        for mo in re.finditer(tok_regex, text):
            kind = mo.lastgroup
            value = mo.group(kind)
            start = mo.start()
            end = mo.end()

            if kind == 'NEWLINE':
                line_start = mo.end()
                line_num += 1
            elif kind == 'MISC':
                pass
            elif kind == 'UNDEFINED':
                raise RuntimeError(f'{value!r} unexpected on line {line_num}')
            else:
                value = value.strip()
                if value == '' or value == None:
                    pass
                else:
                    column = mo.start() - line_start
                    yield self.Token(kind, value, line_num, start, end)
