import re
from strings import *

# Find replication header
# Read replication until next header
# Determine column number
# Determine hardening plan
# Add to hashmap if not present, otherwise increment replication index and assign
# Count initial tripped lines
# Count total tripped lines
# Sum load shed
# Convert to csv and write to file
# Repeat

COLUMN_NUMBER = 4 # Optional explicit definition of hardened line column.
class Analyzer:


    def __init__(self):
        self.hardening_plans = {}
        self.rep_headers_list = []
        self.replication_indicies = []

    def parse(self):
        pass

    def parse_header(self, searchstring):
        myregex = r"=+" + searchstring + r"=+"
        pattern = myregex.compile(myregex)
        return pattern

    def parse_section(self):
        pass

    def search_headers(self, target, string):
        myregex = r"=+" + target + r"=+"
        pattern = re.compile(myregex)
        return pattern.finditer(string)

    def determine_plan(self, replication_text, colnum):
        ln = 0
        for line in replication_text.splitlines():
            print(str(ln) + ": " + line)
            ln += 1


    def build_rep_header_indices(self, log, header):
        headers = self.search_headers(header, log)
        #print(headers)
        #print(log)
        for header in headers:
            self.rep_headers_list.append((header.start(), header.end()))

        print(self.rep_headers_list)
            #print("Start: " + str(header.start()))
            #print("End: " + str(header.end()))
            #print (replication)

    def build_replication_indices(self):
        start = 0
        end = 1
        length = len(self.rep_headers_list)

        for index, header in enumerate(self.rep_headers_list):
            if (index+1 < length):
                rep_index = (header[end]+1, self.rep_headers_list[index+1][start]-1)
                #rep_length = self.rep_headers_list[index+1][start] - header[end]
                self.replication_indicies.append(rep_index)
                print(rep_index)
            #print(header[1])
            #print(self.rep_headers_list[index+1][0])

    def analyze_replication(self, replication_text):
        self.determine_plan(replication_text, COLUMN_NUMBER-1)
        pass

    def read_file(self, filename):
        """
        Read line until we encounter LINE_FAILURE header. Load text until we encounter next header.
        """
        with open(filename, "r") as f:
            log = f.read()
            log_string = str(log)
            self.build_rep_header_indices(log_string, REP_HEADER)
            self.build_replication_indices()
            replication_text = log_string[self.replication_indicies[0][0] : self.replication_indicies[0][1]]
            self.analyze_replication(replication_text)
            #replication = log_string[head_start : head_end]
            #print(str(line))
            #print(LINE_FAILURE)



            """
            replication =
            if header != None:
                #print(line)
                #print(header)
                print(header.string + " " + str(header.endpos))
                #print(line[10:87])
            else:
                print("Fail")

            self.determine_plan(4)
            """

    def aggregate(self):
        pass

    def calculate(self):
        pass

    def display(self):
        pass


analyzer = Analyzer()
analyzer.read_file("cas_results_5.txt")
#print(LINE_FAILURE)
