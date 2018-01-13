import re
from strings import *

class Analyzer:
    hardening_plans = { }
    #def __init__(self)

    def __init__(self):
        self.hardening_plans = { }

    def parse(self):
        pass

    def parse_header(self, searchstring):
        myregex = r"=+" + searchstring + r"=+"
        pattern = re.compile(myregex)
        return pattern

    def parse_section(self):
        pass

    def search_header(self, target, string):
        myregex = r"=+" + target + r"=+"
        pattern = re.compile(myregex)
        return re.search(pattern, string)

    def determine_plan(self, colnum):
        self.hardening_plans[(1, 3, 17)] = 0
        self.hardening_plans[(1)] = 3
        print(self.hardening_plans[(1, 3, 17)])
        print(self.hardening_plans[1])

    def aggregate(self):
        pass

    def calculate(self):
        pass

    def display(self):
        pass

    def read(self, filename):
        """
        Read line until we encounter LINE_FAILURE header. Load text until we encounter next header.
        """
        with open(filename, "r") as f:
            line = f.readline()
            #print(str(line))
            #print(LINE_FAILURE)
            header = self.search_header(LINE_FAILURE, line)
            #header = None
            if header != None:
                #print(line)
                #print(header)
                print(header.string + " " + str(header.endpos))
                #print(line[10:87])
            else:
                print("Fail")

            self.determine_plan(4)


analyzer = Analyzer()
analyzer.read("cas_results_5.txt")
#print(LINE_FAILURE)
