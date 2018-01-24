from processor import *
from tokenizer import *

class Analyzer:

    def __init__(self):
        self.processor = Processor()
        self.tokenizer = Tokenizer()

    def read(self, filename):
        with open(filename, "r") as f:
            log = f.read()
            return str(log)

    def run(self, text, outfile, flag):
        token_collection = self.tokenizer.tokenize(text)
        self.processor.process(token_collection, outfile, flag)
